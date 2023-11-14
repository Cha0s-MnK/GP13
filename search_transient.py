###################
# IMPORT PACKAGES #
###################

import glob
import grand.dataio.root_trees as rt # GRANDLIB
import numpy as np
np.set_printoptions(threshold=np.inf)
import os
import time as wall_time
start_time = wall_time.perf_counter()

#############################################
# SEARCH FOR TRANSIENTS IN GRAND DATA FILES #
#############################################

# GRANDLIB BRANCH: 'dev_io_root'
# This script reads out a GRAND data file in ROOT format. It searches for transient pulses in the time traces of the data 
# file, where each DU is treated separately. Here, a transient is defined as a pulse that exceeds 5 times the STD of the 
# trace, called 5 sigma from this point onwards. For each trace, the number of 5 sigma transients is computed, as well as 
# the number of times this 5 sigma threshold is exceeded during the transient. The maximum peak position of the transient is
# also stored. All information is saved in an NPZ file for each DU.

####################
# DEFINE ARGUMENTS #
####################

channels            = ['X', 'Y', 'Z'] # make dictionaries for easier indexing
data_dir            = 'data/' # path of data directory containing ROOT files to analyze
specific_date       = '20231013/' # specify ROOT files of which date to analyze
num_crossings       = 2 # least number of threshold crossings in a time window
num_threshold       = 5 # trigger threshold of the transient (times of noises)
noises              = [20.0, 20.0, 40.0] # average noise level for 3 ADC channels (ADC counts)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)

# input dates manually
date_list = ['20231011/', '20231012/', '20231013/', '20231014/', '20231015/', '20231016/', '20231017/', '20231018/', 
             '20231019/', '20231020/', '20231027/', '20231028/', '20231029/', '20231030/', '20231031/']

# input DUs manually
du_list = ['1010', '1013', '1016', '1017', '1019', '1020', '1021', '1029', '1031', '1032', '1033', '1035', '1041']

##############################
# GET ROOT FILES AND DU LIST #
##############################

# get ROOT files and their information to analyze
file_list = sorted(glob.glob(os.path.join(data_dir, specific_date, '*.root')))
num_files = len(file_list)

# create an empty set to store unique DUs
du_set = set()

# loop through all files to get used DUs
for file in file_list:
    tadc = rt.TADC(file) # initiate TADC tree of this file
    tadc.get_entry(0) # get the entry from this file
    du_set.update(tadc.get_list_of_all_used_dus()) # update the set with all used DUs from this file

# convert the set to a list in order
du_list = sorted(list(du_set))

# print the list of all used DUs
print(f'\nROOT files contain data from following DUs: \n{du_list}\n')
'''
# get some useful parameters from the data directory
data_dir_split = data_dir.split('/')
array = path_to_data_dir_split[4]
month = path_to_data_dir_split[5]
mode  = path_to_data_dir_split[6]
'''
####################################
# CORE FUNCTIONS TO SEARCH WINDOWS #
####################################

def search_windows(trace, threshold, standard_separation):
    # stop if there are no transients
    exceed_threshold = np.abs(trace) > threshold
    if np.sum(exceed_threshold) < num_crossings:
        return []
    
    # find trace positions where threshold is exceeded
    crossing_ids = np.flatnonzero(exceed_threshold)
    
    # find the separations between consecutive threshold crossings
    crossing_separations = np.diff(crossing_ids)

    # locate pulse indices in threshold crossing indices
    pulse_ids = np.flatnonzero(crossing_separations > standard_separation)
    pulse_ids = np.concatenate(([-1], pulse_ids, [len(crossing_ids)-1]))
    
    # preallocate the return list for time windows
    window_list = []

    # search all transients/pulses
    half_separation = standard_separation // 2
    for i in range(len(pulse_ids)-1):
        # get the start index of current pulse
        start_id = crossing_ids[pulse_ids[i]+1] - half_separation
        start_id = max(0, start_id) # fix the 1st pulse

        # get the stop index of current pulse
        stop_id = crossing_ids[pulse_ids[i+1]] + half_separation
        stop_id = min(len(trace)-1, stop_id) # fix the last pulse

        # check if this time window contains enough crossings
        if np.sum(exceed_threshold[start_id:stop_id+1]) >= num_crossings:
            window_list.append([start_id, stop_id])

    return window_list

def iterative_search(trace, num_threshold, standard_separation):
    # initialize variables
    window_list = []
    cur_window_list = []
    noise_trace = np.copy(trace)

    while True:
        # calculate current threshold using the standard deviation of current noise trace
        threshold = num_threshold * np.std(noise_trace)

        # search for time windows using current threshold
        cur_window_list = search_windows(trace, threshold, standard_separation)

        # check if the number of time windows found is the same as that in the previous iteration
        if len(cur_window_list) == len(window_list):
            return cur_window_list
        else:
            window_list = cur_window_list

            # exclude the windowed parts from the trace to get next noise trace 
            noise_trace = np.copy(trace)
            for start_id, stop_id in cur_window_list:
                noise_trace[start_id:stop_id+1] = 0
            noise_trace = noise_trace[noise_trace != 0]

##################################
# SEARCH TRANSIENTS IN ALL FILES #
##################################

# mask channels to disable channel 0 for GP13 data
channel_mask = np.array([False, True, True, True])

# make dictionaries to avoid looping through all files more than once
windows_list = {channel: {} for channel in channels}
gps_times_list = {}

# loop through all DUs to initiate dictionaries with empty lists
for du in du_list:
    for channel in channels:
        windows_list[channel][du] = []
    gps_times_list[du] = []

# loop through all files
for file_id, file in enumerate(file_list):
    # get info of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()

    # loop over all entries in this file (1 event corresponds to several entries)
    print(f'\n{file_id+1}/{num_files}: Looping over {num_entries} entries in {file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
    
        # get traces and search transients
        traces = np.array(tadc.trace_ch[0])[channel_mask]
        for channel_id, channel in enumerate(channels):
            windows_list[channel][du].append(search_windows(traces[channel_id], num_threshold*np.std(traces[channel_id]), standard_separation))
            #windows_list[channel][du].append(search_windows(traces[channel_id], num_threshold*noises[channel_id], standard_separation))
        gps_times_list[du].append(trawv.gps_time[0])

# loop through all DUs
for du in du_list:
    # save all info into a NPZ file
    npz_dir  = 'result/'
    npz_file = 'DU{}_threshold{}_separation{}_date{}.npz'.format(du, num_threshold, standard_separation, specific_date[:8])
    np.savez(os.path.join(npz_dir, specific_date, npz_file),
             du = du,
             **{f'window{channel}': windows_list[channel][du] for channel in channels},
             gps_times = gps_times_list[du])
    print(f'Saved: {npz_dir}{specific_date}{npz_file}')

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")