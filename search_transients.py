# system
import os
import glob
import argparse
import datetime

# scipy
import numpy as np

# grandlib
import grand.dataio.root_trees as rt

# time
import time
start_time = time.perf_counter()

################################################
# SEARCHING FOR TRANSIENTS IN GRAND DATA FILES #
################################################

# GRANDLIB BRANCH: 'dev_io_root'
# This script reads out a GRAND data file in ROOT format. It searches for transient pulses in the time traces of the data 
# file, where each DU is treated separately. Here, a transient is defined as a pulse that exceeds 5 times the STD of the 
# trace, called 5 sigma from this point onwards. For each trace, the number of 5 sigma transients is computed, as well as 
# the number of times this 5 sigma threshold is exceeded during the transient. The maximum peak position of the transient is
# also stored. All information is saved in an NPZ file for each DU.

###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Search for transient in data traces and save information in NPZ files per DU.")

parser.add_argument('--path_to_data_dir',
                    dest='path_to_data_dir',
                    default='gp13_data/',
                    type=str,
                    help='Specifies the path of the directory containing the ROOT files to analyse.')

parser.add_argument("--file_specific",
                    dest="file_specific",
                    default='*.root',
                    type=str,
                    help="Specify which ROOT files to analyze.")

parser.add_argument('--empty_channel',
                    dest='empty_channel',
                    default=0,
                    type=int,
                    help='Specifies which of the 4 ADC channels is not branched.')

parser.add_argument('--num_threshold',
                    dest='num_threshold',
                    default=5,
                    type=int,
                    help='Specifies the trigger threshold in units of sigma/STDs.')

parser.add_argument('--standard_separation',
                    dest='standard_separation',
                    default=100,
                    type=int,
                    help='Specifies the required separation between two pulses in a trace in units of sample numbers.')                    

parse_args = parser.parse_args()

path_to_data_dir    = parse_args.path_to_data_dir
file_specific       = parse_args.file_specific
empty_channel       = parse_args.empty_channel
num_threshold       = parse_args.num_threshold
standard_separation = parse_args.standard_separation

#####################################
# GET ROOT FILES AND DU INFORMATION #
#####################################

# check if the data directory or the file list is empty
if not os.path.exists(path_to_data_dir):
    raise FileNotFoundError(f'Data directory not found: {path_to_data_dir}')
root_files = sorted(glob.glob(os.path.join(path_to_data_dir, file_specific)))
if not root_files:
    raise ValueError("The provided file list is empty. Please ensure that you provide a non-empty list of ROOT files.")

# initiate TADC tree of the file list to get basic info
total_entries   = 0 # total number of entries across all files
max_dus         = 0 # maximum used DUs
id_max_dus_file = 0 # index of file that has maximum used DUs
for i, file in enumerate(root_files):
    tadc = rt.TADC(file)
    tadc.get_entry(0) # get the entry from current file
    total_entries += tadc.get_number_of_entries()
    cur_dus = len(tadc.get_list_of_all_used_dus()) # get used DUs from current file
    if cur_dus > max_dus:
        max_dus         = cur_dus
        id_max_dus_file = i

# get list of used DUs
tadc = rt.TADC(root_files[id_max_dus_file])
du_list = tadc.get_list_of_all_used_dus()
print('\nFiles contain data from following DUs: {}'.format(du_list))
'''
# Sometimes all four channels are enabled, but one is not branched. For GP13 it is always channel 0
mask_channel = np.array(tadc.adc_enabled_channels_ch[0], dtype=bool)
if np.all(mask_channel):
    mask_channel[empty_channel] = False

# Sometimes there seems to be a bug where there are less than 3 enabled ADC channels.
if len(mask_channel[mask_channel]) < 3:
    mask_channel = np.ones(4, dtype=bool)
    mask_channel[empty_channel] = False
'''
# create a channel mask to enable 3 ADC channels and disable channel 0 for GP13
mask_channel = np.ones(4, dtype=bool)
mask_channel[empty_channel] = False
'''
# get some useful parameters from the data directory
path_to_data_dir_split = path_to_data_dir.split('/')
array = path_to_data_dir_split[4]
month = path_to_data_dir_split[5]
mode  = path_to_data_dir_split[6]

# define the NPZ directory and file name template to store the computed PSDs
npz_dir           = 'data/{}/{}/{}/transient_search/'.format(array, month, mode)
npz_file_template = 'transient_search_{}_{}_{}_du_{}_thresh_{}_separ_{}.npz'.format(array, month, mode, '{}', thresh, standard_separation)
'''
# define the NPZ directory and file name template to store the computed PSDs
npz_dir           = 'result/'
npz_file_template = 'du_{}_threshold_{}_separation_{}.npz'.format('{}', num_threshold, standard_separation)

def search_transients(trace,
                      num_threshold=5,
                      standard_separation=100): 
    '''
    Identify transients of a signal where the signal exceeds a certain threshold.

    Parameters:
    trace (array-like): the signal to be analyzed
    num_threshold (float): a multiplier for the STD to determine the threshold for identifying transients.
    standard_separation (int): a minimum separation between transients in units of samples.

    Returns:
    windows: a list of windows, each containing the start and end points of a transient.
    crossings: a list of crossings, each containing the values of the signal during a transient.
    '''
    # create the return lists
    windows   = [] # time windows
    crossings = [] # threshold crossings

    # find trace positions where N sigma/STD is exceeded
    threshold        = num_threshold * np.std(trace)
    exceed_threshold = np.abs(trace) > threshold
    ids_crossings    = np.nonzero(exceed_threshold)[0]

    # stop if there are no transients
    if not np.any(exceed_threshold):
        return windows, crossings

    # find the separations between consecutive threshold crossings
    separations_crossings = np.diff(ids_crossings)

    # locate where two transients/pulses are separated
    ids_pulses = np.concatenate((np.nonzero(separations_crossings > standard_separation)[0], [len(ids_crossings)-1]))
    
    # search all transients/pulses
    half_separation = standard_separation / 2
    for i, id in enumerate(ids_pulses):
        # get the beginning of current pulse
        if i == 0:
            begin_pulse = np.max([0, ids_crossings[0] - half_separation])
        else:
            begin_pulse = end_pulse - standard_separation + separations_crossings[ids_pulses[i - 1]]

        # get the end of current pulse
        end_pulse = np.min([len(trace)-1, ids_crossings[id] + half_separation])

        # count how many crossings (in units of sigma/STD) are in the time window
        abs_trace_window = np.abs(trace[int(begin_pulse) : int(end_pulse)+1])
        crossing         = [signal * num_threshold / threshold for signal in abs_trace_window if signal > threshold]

        # add info of current pulse to the return lists
        windows.append([int(begin_pulse),int(end_pulse)])
        crossings.append(list(crossing))

    return windows, crossings

#########################################
# COMPUTE AND PLOT MEAN PSD FOR EACH DU #
#########################################

# make dictionaries so that you don't have to loop over all files each time you want to analyze a different DU
channels  = ['x', 'y', 'z']
traces = {channel: {} for channel in channels}
windows = {channel: {} for channel in channels}
crossings = {channel: {} for channel in channels}
gps_time = {}
for du in du_list:
    for channel in channels:
        windows[channel][du]   = np.zeros((total_entries), dtype=object) # save time windows of traces (unit = sample number)
        crossings[channel][du] = np.zeros((total_entries), dtype=object) # save threshold crossings of transients (unit = sigma/STD)
    gps_time[du] = np.zeros((total_entries), dtype=int) # save GPS times

# loop over all files in run
num_files = len(root_files)
id_entry  = 0
for i, root_file in enumerate(root_files[:5]):
    tadc  = rt.TADC(root_file)
    trawv = rt.TRawVoltage(root_file)
    num_entries = tadc.get_number_of_entries()

    # loop over all events in file
    print('\n{}/{}: Looping over {} events in {}'.format(i+1, num_files, num_entries, root_file))
    for entry in range(num_entries):
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # use the DU as key for the dictionaries
        # Assume that there are no coincident events, i.e. only one DU (and trace etc) per event.
        du = tadc.du_id[0]
    
        # get traces
        new_traces = np.array(tadc.trace_ch[0], dtype=object)[mask_channel]
        for channel, new_trace in zip(channels, new_traces):
            traces[channel] = new_trace
        '''
        # check if sample numbers are enough
        if any(len(channel) != num_samples for channel in new_traces[:3]):
            print('\nWARNING for entry {}: {} samples in trace != (3,{})'.format(entry, new_traces.shape, num_samples) + '\nSkipping...')
            id_entry += 1
            continue
        
        # check if data is between 10-55ºC
        if trawv.gps_temp[0] > 55 or trawv.gps_temp[0] < 10:
            print('\nWARNING for entry {}: gps_temp = {} not between 10 and 55 ºC'.format(entry,trawv.gps_temp[0]) + '\nSkipping...')
            id_entry += 1
            continue
        '''
        for channel in channels:
            #'''
            new_windows, new_crossings = search_transients(traces[channel],
                                                           num_threshold=num_threshold,
                                                           standard_separation=standard_separation)
            '''
            new_windows, new_crossings = search_transients_origin(traces[channel],
                                                                  thresh=num_threshold,
                                                                  separation_pulses=standard_separation)
            '''
            windows[channel][du][id_entry] = new_windows
            crossings[channel][du][id_entry] = new_crossings
        gps_time[du][id_entry] = trawv.gps_time[0]

        id_entry += 1

for du in du_list:
    # get rid of events that are not filled because of bad quality data
    mask_filled = np.where(windows['x'][du] != 0)
    for channel in channels:
        windows[channel][du] = windows[channel][du][mask_filled]
        crossings[channel][du] = crossings[channel][du][mask_filled]
    gps_time[du] = gps_time[du][mask_filled]

    # save total power in NPZ file
    npz_file = npz_file_template.format(du)
    np.savez(f'{npz_dir}{npz_file}',
             du = du,
             window_x = windows['x'][du],
             window_y = windows['y'][du],
             window_z = windows['z'][du],
             crossing_x = crossings['x'][du],
             crossing_y = crossings['y'][du],
             crossing_z = crossings['z'][du],
             gps_time = gps_time[du])
    print(f'\nSaved NPZ file: {npz_dir}{npz_file}')

# record the running time
end_time = time.perf_counter()
run_time = end_time - start_time
print(f"Whole program executed in: {run_time} seconds")

# read the selected NPZ file
npz_file = np.load('result/du_1013_threshold_5_separation_100.npz', allow_pickle=True)

# use a list comprehension to filter out empty lists
windows_x = np.array([a_list for a_list in npz_file['window_x'] if len(a_list) > 0])
print(windows_x)