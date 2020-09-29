#!/usr/bin/env python3

from os import mkdir, remove
from os.path import isdir


import numpy as np

from LoLIM.utilities import processed_data_dir
from LoLIM.IO.raw_tbb_IO import filePaths_by_stationName, MultiFile_Dal1, read_station_delays, read_bad_antennas, read_antenna_delays
from LoLIM.IO.metadata import ITRF_to_geoditic
from LoLIM.findRFI import window_and_filter
from LoLIM.interferometry import read_interferometric_PSE as R_IPSE
from LoLIM.getTrace_fromLoc import getTrace_fromLoc
from LoLIM.signal_processing import remove_saturation, upsample_and_correlate, parabolic_fit, num_double_zeros, data_cut_inspan, locate_data_loss
from LoLIM.signal_processing import parabolic_fitter

import sys 
sys.path.append('..')
            
from LMA_window_data import writer_manager, antenna_symbols, window_data

if __name__ == "__main__":
    from LoLIM import utilities
    utilities.default_raw_data_loc = "/home/brian/KAP_data_link/lightning_data"
    utilities.default_processed_data_loc = "/home/brian/processed_files"
    
    stations = [('CS002',5), ('RS205',4),  ('RS306',4), ('RS406',5), ('RS106',5), ('RS407',5), ('RS409',5), ('RS208',5), ('RS508',5), ('RS310',5)]  
    
    
    
    start_sample = 500*(2**16)
    end_sample = 6000*(2**16)
    window_length_time = 80E-6
    amp_tresh = 10
    blocksize = 2**16
    
    timeID = 'D20190424T194432.504Z'
    polarization_flips = 'polarization_flips.txt'
    bad_antennas = 'bad_antennas.txt'
    additional_antenna_delays = 'ant_delays.txt'
    station_delays_fname = 'station_delays.txt'
    pol_flips_are_bad = True
    
    num_dataLoss_zeros = 10
    
    hann_window_fraction = 0.1
    lower_frequency = 60.0E6
    upper_frequency = 66.0E6
    
    
    ### sanity check
    num_used = np.sum([ d[1] for d in  stations])
    if num_used > len(antenna_symbols):
        print("ERROR: too many antennas requested!")
        quit()
    
    
    processed_data_folder = processed_data_dir( timeID )
    station_timing_offsets = read_station_delays( processed_data_folder+'/'+ station_delays_fname )
    raw_fpaths = filePaths_by_stationName( timeID )
    
    
    data_filter = window_and_filter(blocksize=blocksize, half_window_percent=hann_window_fraction, 
                                    lower_filter=lower_frequency, upper_filter=upper_frequency)
    
    writer_manager = writer_manager('/home/brian/paper_scripts/LMA_comparison/2019_Traditional80us/text_out')
    for sname,num_ants in stations:
        print(sname)
        TBB_data = MultiFile_Dal1(raw_fpaths[sname], polarization_flips = processed_data_folder+'/'+polarization_flips, 
          bad_antennas = processed_data_folder+'/'+bad_antennas, additional_ant_delays = processed_data_folder+'/'+additional_antenna_delays,
          pol_flips_are_bad = pol_flips_are_bad)
                
        TBB_data.set_station_delay( station_timing_offsets[sname] )
        TBB_data.find_and_set_polarization_delay()
        
        edge_length = int(hann_window_fraction*blocksize) + 10
        
#        (TBB_data, filter, edge_length, window_length_time, start_sample_number, end_sample_number, amp_tresh, num_dataLoss_zeros, max_num_antennas, data_writer_manager)
        f = window_data( TBB_data, data_filter, edge_length, window_length_time, start_sample, end_sample, amp_tresh, num_dataLoss_zeros, num_ants, writer_manager, 
                        inject_T_noise=None)
        print("  fraction without dataloss:",f)
        
    writer_manager.write_loc_file()
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    