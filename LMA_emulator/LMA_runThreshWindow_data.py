#!/usr/bin/env python3

from os import mkdir, remove
from os.path import isdir
from datetime import datetime, timezone, timedelta

#import struct

import numpy as np

from scipy.special import erf
from scipy.optimize import root_scalar

from LoLIM.utilities import processed_data_dir, v_air
from LoLIM.IO.raw_tbb_IO import filePaths_by_stationName, MultiFile_Dal1, read_station_delays, read_antenna_pol_flips, read_bad_antennas, read_antenna_delays
from LoLIM.IO.metadata import ITRF_to_geoditic
from LoLIM.findRFI import window_and_filter
from LoLIM.interferometry import read_interferometric_PSE as R_IPSE
from LoLIM.getTrace_fromLoc import getTrace_fromLoc
from LoLIM.signal_processing import remove_saturation, upsample_and_correlate, parabolic_fit, num_double_zeros, data_cut_inspan, locate_data_loss
from LoLIM.signal_processing import parabolic_fitter

from LMA_window_data import writer_manager, antenna_symbols




def window_data(TBB_data, filter, edge_length,  start_sample_number, end_sample_number, amp_tresh, 
        num_dataLoss_zeros, max_num_antennas, data_writer_manager, inject_T_noise=None, 
        bin_width_time=10e-6, rate=1.0/80E-6, num_bins=40, shift_precent=0.1):
    
    if inject_T_noise is not None:
        clock_noise = np.random.normal(scale=inject_T_noise)
    else:
        clock_noise = 0.0
        
    ITRFantenna_locations = TBB_data.get_ITRF_antenna_positions()
    antenna_names = TBB_data.get_antenna_names()
    sname = TBB_data.get_station_name()
    posix_timestamp = TBB_data.get_timestamp()
    
    num_ants = len(ITRFantenna_locations)
    num_station_antennas = len(ITRFantenna_locations)
    num_ants = min(int(num_station_antennas / 2), max_num_antennas)
    antennas_to_use = np.arange(num_ants) * 2
    
    writers = []
    for antI in antennas_to_use:
        name = antenna_names[antI]
        lat_lon_alt = ITRF_to_geoditic(ITRFantenna_locations[antI])
        writers.append(data_writer_manager.get_next_writer(name, lat_lon_alt, posix_timestamp))

    antenna_start_times = TBB_data.get_time_from_second()
    antenna_start_times += clock_noise
    
    approx_start_T = np.average(antenna_start_times) + start_sample_number * 5e-09
    approx_end_T = np.average(antenna_start_times) + end_sample_number * 5e-09
    
    window_index = int(approx_start_T / bin_width_time)
    blocksize = filter.blocksize
    data_loss_segments = [[]] * num_station_antennas
    workspace = np.empty(blocksize, dtype=(np.complex))
    
    data_blocks = np.empty((num_station_antennas, blocksize), dtype=(np.double))
    current_samples = np.empty(num_station_antennas, dtype=(np.int))

    def load_data_block(ant_i, start_T):
        sample_number = int((start_T - antenna_start_times[ant_i]) / 5e-09)
        sample_number -= edge_length
        TMP = TBB_data.get_data(sample_number, blocksize, antenna_index=ant_i)
        data_loss_spans, DL = locate_data_loss(TMP, num_dataLoss_zeros)
        workspace[:] = TMP
        workspace[:] = filter.filter(workspace)
        np.abs(workspace, out=(data_blocks[ant_i, :]))
        data_loss_segments[ant_i] = data_loss_spans
        current_samples[ant_i] = sample_number

    for ant_i in antennas_to_use:
        load_data_block(ant_i, window_index * bin_width_time)

    print_width = 100
    windows_total = 0
    windows_found = 0
    
    current_threshold = amp_tresh
    current_num_pulses = 0
    current_num_bins = 0
    previous_rate = None
    
    last_pulse_window_index = None
    previous_peak_time = None
    window_end_T = 0
    while window_end_T < approx_end_T:
        windows_total += 1
        window_start_T = window_index * bin_width_time
        window_end_T = window_start_T + bin_width_time
        
        if not window_index % print_width:
            print(sname, window_index, (window_start_T - approx_start_T) / (approx_end_T - approx_start_T))
            if previous_rate is not None:
                print("  ", current_threshold, previous_rate/rate)
           

        for ant_i, writer in zip(antennas_to_use, writers):
            writer.write_window_boundry(window_start_T)
            
        
            local_start_index = int((window_start_T - antenna_start_times[ant_i]) / 5e-09) - current_samples[ant_i]
            actual_start_time = (local_start_index + current_samples[ant_i]) * 5e-09 + antenna_start_times[ant_i]
            actual_index = int(actual_start_time / bin_width_time)
            while actual_index != window_index:
                local_start_index -= np.sign(actual_index - window_index)
                actual_start_time = (local_start_index + current_samples[ant_i]) * 5e-09 + antenna_start_times[ant_i]
                actual_index = int(actual_start_time / bin_width_time)


            local_end_index = local_start_index + int(bin_width_time / 5e-09)
            if local_end_index > blocksize - edge_length:
                load_data_block(ant_i, window_start_T)
                local_start_index = int((window_start_T - antenna_start_times[ant_i]) / 5e-09) - current_samples[ant_i]
                local_end_index = local_start_index + int(bin_width_time / 5e-09)
                
            has_dataloss = data_cut_inspan(data_loss_segments[ant_i], local_start_index, local_end_index)
            
            if not has_dataloss:
                windows_found += 1
                window = data_blocks[ant_i, local_start_index + 1:local_end_index]
                max_point = np.argmax(window)
                peak_time = (max_point + 1 + local_start_index + current_samples[ant_i]) * 5e-09 + antenna_start_times[ant_i]
                amp = window[max_point]
                
                if amp > current_threshold:
                    writer.write_pulse(peak_time)
                    current_num_pulses += 1
                    
                    pulse_index = int((peak_time - int(peak_time)) / bin_width_time)
                    if last_pulse_window_index is not None:
                        if pulse_index == last_pulse_window_index:
                            print('windowing error', previous_peak_time, peak_time)
                            print(actual_start_time, max_point)
                    last_pulse_window_index = pulse_index
                    previous_peak_time = peak_time
                    
                ## update our running variables
                current_num_bins += 1
                
                if current_num_bins == num_bins:
                    current_rate = current_num_pulses/(current_num_bins*bin_width_time)
                    previous_rate = current_rate
                    current_num_bins = 0
                    current_num_pulses = 0
                    
                    if current_rate > rate:
                        current_threshold *= (1.0+shift_precent)
                            
                    elif current_rate < rate:
                        current_threshold *= (1.0-shift_precent)
                        if current_threshold<amp_tresh:
                            current_threshold = amp_tresh
                        
                    
                break
                    
                    

        window_index += 1

    for w in writers:
        w.finalize()

    return windows_found / windows_total



def choose_S(f):
    inv_sqrt2 = 1.0/np.sqrt(2)
    def func(S):
        return erf(S*inv_sqrt2) - f
    
    res = root_scalar( func, x0=1, method='brentq', bracket=(0,10) )
    print( res )
    
if __name__=="__main__":
    choose_S( 7.0/8.0 )
    

    
    

    
    
    
    