#!/usr/bin/env python3

# uncompyle6 version 3.6.4
# Python bytecode 3.6 (3379)
# Decompiled from: Python 3.6.8 |Anaconda, Inc.| (default, Dec 30 2018, 01:22:34) 
# [GCC 7.3.0]
# Embedded file name: ../LMA_window_data.py
# Compiled at: 2020-02-26 11:28:13
# Size of source mod 2**32: 20546 bytes
from os import mkdir, remove
from os.path import isdir
from datetime import datetime, timezone, timedelta
import numpy as np
from LoLIM.utilities import processed_data_dir, v_air
from LoLIM.IO.raw_tbb_IO import filePaths_by_stationName, MultiFile_Dal1, read_station_delays, read_antenna_pol_flips, read_bad_antennas, read_antenna_delays
from LoLIM.IO.metadata import ITRF_to_geoditic
from LoLIM.findRFI import window_and_filter
from LoLIM.interferometry import read_interferometric_PSE as R_IPSE
from LoLIM.getTrace_fromLoc import getTrace_fromLoc
from LoLIM.signal_processing import remove_saturation, upsample_and_correlate, parabolic_fit, num_double_zeros, data_cut_inspan, locate_data_loss
from LoLIM.signal_processing import parabolic_fitter
antenna_symbols = [
 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'x', 'y', 'z']

class writer_manager:

    class antenna_writer:

        def __init__(self, file_name, posix_timestamp, ant_name):
            self.file_name = file_name
            self.ant_name = ant_name
            self.fout = open(self.file_name, 'w')
            self.window_file_name = self.file_name.replace('.txt', '.win')
            self.window_fout = open(self.window_file_name, 'w')
            self.timestamp = datetime.fromtimestamp(posix_timestamp, tz=(timezone(timedelta(seconds=0))))
            self.num_pulses = 0

        def write_pulse(self, time_from_posix_stamp, thresh=50, amp=1):
            seconds = int(time_from_posix_stamp)
            subseconds = time_from_posix_stamp - seconds
            TS = self.timestamp + timedelta(seconds=seconds)
            self.fout.write(TS.strftime('%Y-%m-%dT%H:%M:%S.'))
            self.fout.write('%09d' % int(subseconds * 1000000000))
            self.fout.write(' ')
            self.fout.write(str(int(thresh)))
            self.fout.write(' ')
            self.fout.write(str(int(amp)))
            self.fout.write('\n')
            self.num_pulses += 1

        def write_window_boundry(self, window_from_posix_stamp):
            seconds = int(window_from_posix_stamp)
            subseconds = window_from_posix_stamp - seconds
            TS = self.timestamp + timedelta(seconds=seconds)
            self.window_fout.write(TS.strftime('%Y-%m-%dT%H:%M:%S.'))
            self.window_fout.write('%09d' % int(subseconds * 1000000000))
            self.window_fout.write('\n')

        def finalize(self):
            self.fout.close()
            if self.num_pulses == 0:
                remove(self.file_name)

    def __init__(self, folder, save_timing_errors=False):
        self.folder = folder
        self.next_ant_i = 0
        self.network_name = 'LOFAR_LMA_EMULATOR'
        self.current_group = 1
        self.antenna_names = []
        self.antenna_IDs = []
        self.antennas_positions = []
        self.antenna_timing_errors = []
        
        self.save_timing_errors = save_timing_errors
        if save_timing_errors:
            self.timing_errors_out = open(self.folder + '/timing_error.err', 'w')

    def get_next_writer(self, ant_name, lat_lon_alt, posix_timeID, time_error=0.0):
        """call this when starting to save new antenna"""
        if self.next_ant_i >= len(antenna_symbols):
            print('can only have', len(antenna_symbols), 'antennas.')
            quit()
        antenna_letter = antenna_symbols[self.next_ant_i]
        self.next_ant_i += 1
        new_writer = self.antenna_writer(self.folder + '/' + antenna_letter + '.txt', posix_timeID, ant_name)
        new_writer.fout.write(antenna_letter)
        new_writer.fout.write('\n')
        new_writer.fout.write(self.network_name)
        new_writer.fout.write('\n')
        new_writer.fout.write(ant_name)
        new_writer.fout.write('\n')
        new_writer.fout.write(str(self.current_group))
        new_writer.fout.write('\n')
        new_writer.fout.write(str(lat_lon_alt[0]))
        new_writer.fout.write('\n')
        new_writer.fout.write(str(lat_lon_alt[1]))
        new_writer.fout.write('\n')
        new_writer.fout.write(str(lat_lon_alt[2]))
        new_writer.fout.write('\n')
        if self.next_ant_i == len(antenna_symbols):
            self.current_group += 1
            self.next_ant_i = 0
        self.antenna_names.append(ant_name)
        self.antenna_IDs.append(antenna_letter)
        self.antennas_positions.append(lat_lon_alt)
        self.antenna_timing_errors.append( float(time_error) )
        
        if self.save_timing_errors:
            self.timing_errors_out.write( antenna_letter )
            self.timing_errors_out.write( ' ' )
            self.timing_errors_out.write( ant_name )
            self.timing_errors_out.write( ' ' )
            self.timing_errors_out.write( str(time_error) )
            self.timing_errors_out.write( '\n' )
        
        
        return new_writer

    def write_loc_file(self):
        if self.current_group != 1:
            print('ERROR!')
            quit()
        fname = self.folder + '/LOCATIONS.loc'
        fout = open(fname, 'w')
        fout.write('# LOFAR location file\n# data format:#    station mnemonic\n#    station ID\n#    latitude (deg)\n#    longitude (deg)\n#    altitude AMSL (m)\n#    cable delay (ns)\n#    LMA board revision #\n#    LMA receiver channel\n# Network Name\n')
        fout.write(self.network_name)
        fout.write('\n# Center of Array\n')
        fout.write(str(self.antennas_positions[0][0]))
        fout.write('\n')
        fout.write(str(self.antennas_positions[0][1]))
        fout.write('\n')
        fout.write(str(self.antennas_positions[0][2]))
        fout.write('\n')
        for name, ID, loc, terr in zip(self.antenna_names, self.antenna_IDs, self.antennas_positions, self.antenna_timing_errors):
            fout.write('# station\n')
            fout.write(name)
            fout.write('\n')
            fout.write(ID)
            fout.write('\n')
            fout.write(str(loc[0]))
            fout.write('\n')
            fout.write(str(loc[1]))
            fout.write('\n')
            fout.write(str(loc[2]))
            fout.write('\n')
            fout.write( str(  int(-terr*10e10)/10 ) ) ## station delay
            fout.write('\n')
            fout.write('3\n3\n') ##  board, LMA reciever

        fout.write('END')


def window_data(TBB_data, filter, edge_length, window_length_time, start_sample_number, end_sample_number, amp_tresh, num_dataLoss_zeros, max_num_antennas, data_writer_manager, inject_T_noise=1e-08):
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
    num_ants = min(int(num_station_antennas / 2), max_num_antennas)
    antennas_to_use = np.arange(num_ants) * 2
    writers = []
    for antI in antennas_to_use:
        name = antenna_names[antI]
        lat_lon_alt = ITRF_to_geoditic(ITRFantenna_locations[antI])
        writers.append(data_writer_manager.get_next_writer(name, lat_lon_alt, posix_timestamp, clock_noise)) ## this puts clock noise in calibration file

    antenna_start_times = TBB_data.get_time_from_second()
    approx_start_T = np.average(antenna_start_times) + start_sample_number * 5e-09
    approx_end_T = np.average(antenna_start_times) + end_sample_number * 5e-09
    window_index = int(approx_start_T / window_length_time)
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
        load_data_block(ant_i, window_index * window_length_time)

    print_width = 100
    windows_total = 0
    windows_found = 0
    last_pulse_window_index = None
    previous_peak_time = None
    window_end_T = 0
    while window_end_T < approx_end_T:
        windows_total += 1
        window_start_T = window_index * window_length_time
        window_end_T = window_start_T + window_length_time
        if not window_index % print_width:
            print(sname, window_index, (window_start_T - approx_start_T) / (approx_end_T - approx_start_T))
            
            
        for ant_i, writer in zip(antennas_to_use, writers):
            writer.write_window_boundry(window_start_T)
            
            local_start_index = int((window_start_T - antenna_start_times[ant_i]) / 5e-09) - current_samples[ant_i]
            actual_start_time = (local_start_index + current_samples[ant_i]) * 5e-09 + antenna_start_times[ant_i]
            actual_index = int(actual_start_time / window_length_time)
            while actual_index != window_index:
                local_start_index -= np.sign(actual_index - window_index)
                actual_start_time = (local_start_index + current_samples[ant_i]) * 5e-09 + antenna_start_times[ant_i]
                actual_index = int(actual_start_time / window_length_time)

            local_end_index = local_start_index + int(window_length_time / 5e-09)
            if local_end_index > blocksize - edge_length:
                load_data_block(ant_i, window_start_T)
                local_start_index = int((window_start_T - antenna_start_times[ant_i]) / 5e-09) - current_samples[ant_i]
                local_end_index = local_start_index + int(window_length_time / 5e-09)
                
            has_dataloss = data_cut_inspan(data_loss_segments[ant_i], local_start_index, local_end_index)
            
            if not has_dataloss:
                windows_found += 1
                window = data_blocks[ant_i, local_start_index + 1:local_end_index]
                max_point = np.argmax(window)
                peak_time = (max_point + 1 + local_start_index + current_samples[ant_i]) * 5e-09 + antenna_start_times[ant_i]
                amp = window[max_point]
                if amp > amp_tresh:
                    writer.write_pulse( peak_time )
                    pulse_index = int((peak_time - int(peak_time)) / window_length_time)
                    if last_pulse_window_index is not None:
                        if pulse_index == last_pulse_window_index:
                            print('windowing error', previous_peak_time, peak_time)
                            print(actual_start_time, max_point)
                    last_pulse_window_index = pulse_index
                    previous_peak_time = peak_time
                break

        window_index += 1

    for w in writers:
        w.finalize()

    return windows_found / windows_total


# if __name__ == '__main__':
#     from LoLIM import utilities
#     utilities.default_raw_data_loc = '/home/brian/KAP_data_link/lightning_data'
#     utilities.default_processed_data_loc = '/home/brian/processed_files'
#     out_folder = './windowed_normal_80_60_66'
#     stations = [('CS002', 5), ('RS205', 5), ('RS306', 5), ('RS406', 5), ('RS307', 5), ('RS407', 5), ('RS409', 5), ('RS208', 5), ('RS508', 4), ('RS310', 4)]
#     start_sample = 196608000
#     end_sample = 376832000
#     window_length_time = 8e-05
#     amp_tresh = 10
#     blocksize = 65536
#     timeID = 'D20180813T153001.413Z'
#     polarization_flips = 'polarization_flips.txt'
#     bad_antennas = 'bad_antennas.txt'
#     additional_antenna_delays = 'ant_delays.txt'
#     station_delays_fname = 'station_delays_2.txt'
#     pol_flips_are_bad = True
#     num_dataLoss_zeros = 10
#     hann_window_fraction = 0.1
#     lower_frequency = 60000000.0
#     upper_frequency = 66000000.0
#     num_used = np.sum([d[1] for d in stations])
#     if num_used > len(antenna_symbols):
#         print('ERROR: too many antennas requested!')
#         quit()
#     processed_data_folder = processed_data_dir(timeID)
#     station_timing_offsets = read_station_delays(processed_data_folder + '/' + station_delays_fname)
#     raw_fpaths = filePaths_by_stationName(timeID)
#     data_filter = window_and_filter(blocksize=blocksize, half_window_percent=hann_window_fraction, lower_filter=lower_frequency,
#       upper_filter=upper_frequency)
#     writer_manager = writer_manager('/home/brian/paper_scripts/LMA_comparison/default_LMA_emulation/text_out2')
#     for sname, num_ants in stations:
#         print(sname)
#         TBB_data = MultiFile_Dal1((raw_fpaths[sname]), polarization_flips=(processed_data_folder + '/' + polarization_flips), bad_antennas=(processed_data_folder + '/' + bad_antennas),
#           additional_ant_delays=(processed_data_folder + '/' + additional_antenna_delays),
#           pol_flips_are_bad=pol_flips_are_bad)
#         TBB_data.set_station_delay(station_timing_offsets[sname])
#         TBB_data.find_and_set_polarization_delay()
#         edge_length = int(hann_window_fraction * blocksize) + 10
#         f = window_data(TBB_data, data_filter, edge_length, window_length_time, start_sample, end_sample, amp_tresh, num_dataLoss_zeros, num_ants, writer_manager)
#         print('  fraction without dataloss:', f)

#     writer_manager.write_loc_file()
# # okay decompiling LMA_window_data.cpython-36.pyc
