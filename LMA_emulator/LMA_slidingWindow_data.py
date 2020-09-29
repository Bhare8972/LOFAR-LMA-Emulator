#!/usr/bin/env python3

from os import mkdir, remove
from os.path import isdir
from datetime import datetime, timezone, timedelta

#import struct

import numpy as np

from LoLIM.utilities import processed_data_dir, v_air
from LoLIM.IO.raw_tbb_IO import filePaths_by_stationName, MultiFile_Dal1, read_station_delays, read_antenna_pol_flips, read_bad_antennas, read_antenna_delays
from LoLIM.IO.metadata import ITRF_to_geoditic
from LoLIM.findRFI import window_and_filter
from LoLIM.interferometry import read_interferometric_PSE as R_IPSE
from LoLIM.getTrace_fromLoc import getTrace_fromLoc
from LoLIM.signal_processing import remove_saturation, upsample_and_correlate, parabolic_fit, num_double_zeros, data_cut_inspan, locate_data_loss
from LoLIM.signal_processing import parabolic_fitter

from LMA_window_data import writer_manager, antenna_symbols

def window_data(TBB_data, filter, edge_length, start_sample_number, end_sample_number, amp_tresh, num_dataLoss_zeros, max_num_antennas, data_writer_manager, 
                inject_T_noise=None):
    
    ## note, the clock noise is one value per all pulses on station
    if inject_T_noise is not None:
        clock_noise = np.random.normal(scale=inject_T_noise)
    else:
        clock_noise = 0.0
        
    ITRFantenna_locations = TBB_data.get_ITRF_antenna_positions()
    antenna_names = TBB_data.get_antenna_names()
    sname = TBB_data.get_station_name()
    num_ants = len(ITRFantenna_locations)
    posix_timestamp = TBB_data.get_timestamp()
    
    num_station_antennas = len(ITRFantenna_locations)
    num_ants = min(int(num_station_antennas/2),max_num_antennas )
    antennas_to_use = np.arange( num_ants )*2
    
    writers = []
    for antI in antennas_to_use:
        name = antenna_names[ antI ]
        lat_lon_alt = ITRF_to_geoditic( ITRFantenna_locations[antI] )
        writers.append( data_writer_manager.get_next_writer(name, lat_lon_alt, posix_timestamp ) )
            
    antenna_start_times = TBB_data.get_time_from_second()
    antenna_start_times += clock_noise
    
    
    
    #### allocate memory ####
    blocksize = filter.blocksize
    data_loss_segments = [ [] ]*num_station_antennas ## note: length of all antennas in station, not just antennas to load
    workspace = np.empty( blocksize, dtype=np.complex )
    data_blocks = np.empty( (num_station_antennas,blocksize), dtype=np.double )## note: length of all antennas in station, not just antennas to load
    current_samples = np.empty( num_station_antennas, dtype=np.int ) ## note: length of all antennas in station, not just antennas to load
    
    #### initialize data
    def load_data_block(ant_i, sample_number):
        sample_number -= edge_length
        
        TMP = TBB_data.get_data(sample_number, blocksize, antenna_index=ant_i )
        data_loss_spans, DL = locate_data_loss(TMP, num_dataLoss_zeros)
        workspace[:] = TMP
        workspace[:] = filter.filter( workspace )
        np.abs(workspace, out = data_blocks[ant_i,:])
        
        data_loss_segments[ant_i] = data_loss_spans
        current_samples[ant_i] = sample_number
        
        return  data_blocks[ant_i, edge_length:-edge_length], len(data_loss_spans)>0
        
    for ant_i in antennas_to_use:
        load_data_block( ant_i, start_sample_number )
    
    print_width = 10 ## number of blocks
    
    bin_width = int( (80e-6)/(5e-9) )
    half_bin_width = int( bin_width/2 )
    
    nominal_blocksize = blocksize - 2*edge_length - bin_width
    number_blocks = ((end_sample_number-start_sample_number)/nominal_blocksize) + 1
    
    blocks_total = 0
    blocks_found = 0
    
    last_pulse_10us_index = None
    
    current_sample = start_sample_number
    last_peak_time = None
    while current_sample<end_sample_number:
        blocks_total += 1
        
        if not blocks_total % print_width:
            print(sname, blocks_total, blocks_total/number_blocks  )
        
        for ant_i,writer in zip(antennas_to_use,writers):
            
            data, has_dataloss = load_data_block(ant_i, current_sample)
            
            if has_dataloss:
                continue
            
            local_start_time = current_sample*5.0E-9 + antenna_start_times[ant_i]
            blocks_found += 1
            
            start_i = half_bin_width
            if (last_peak_time is not None) and ((local_start_time + start_i*5e-9) < last_peak_time ):
                start_i = int( (last_peak_time-local_start_time)/5e-9 ) + 1
            
            for i in range(half_bin_width, len(data)-half_bin_width):
                if data[i] > amp_tresh:
                
                    window = data[i-half_bin_width:i+half_bin_width]
                    AM = np.argmax( window )
                    if AM == half_bin_width :
                        peak_time = local_start_time + i*5.0E-9
                        
                        pulse_10us_index = int( (peak_time-int(peak_time))/(10e-6) )
                        if (last_pulse_10us_index is None or pulse_10us_index != last_pulse_10us_index) and (last_peak_time is None or (peak_time-last_peak_time)>(40e-6) ) :
                            
                            if (last_peak_time is not None) and peak_time<last_peak_time:
                                print("giggles McDukes!")
                                quit()
                            
                            writer.write_pulse( peak_time )
                            last_peak_time = peak_time
                            
                        last_pulse_10us_index = pulse_10us_index
            break
        current_sample += nominal_blocksize
    
    for w in writers:
        w.finalize()
    
    return blocks_found/blocks_total
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    