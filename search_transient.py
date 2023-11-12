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

# input date list manually
date_list = ['20231011/', '20231012/', '20231013/', '20231014/', '20231015/', '20231016/', '20231017/', '20231018/', 
             '20231019/', '20231020/', '20231027/', '20231028/', '20231029/', '20231030/', '20231031/']

channel_list        = ['X', 'Y', 'Z'] # make dictionaries for easier indexing
data_dir            = 'data/' # path of data directory containing ROOT files to analyze
specific_date       = '20231031/' # specify ROOT files of which date to analyze
num_threshold       = 5 # trigger threshold of the transient (times of noises)
noise_list          = [20.0, 20.0, 40.0] # average noise level for 3 ADC channels (ADC counts)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)

##############################
# GET ROOT FILES AND DU LIST #
##############################

# get ROOT files and their information to analyze
file_list = sorted(glob.glob(os.path.join(data_dir, specific_date, '*.root')))
num_files = len(file_list)

# get list of all used DUs

# create an empty set to store unique DUs
du_set = set()

# loop through all files to get used DUs
for file in file_list:
    tadc = rt.TADC(file) # initiate TADC tree of this file
    tadc.get_entry(0) # get the entry from this file
    used_dus = tadc.get_list_of_all_used_dus() # get used DUs from this file
    du_set.update(used_dus) # update the set

# convert the set to a list in order
du_list = sorted(list(du_set))

# print the list of all used DUs
print(f'\nROOT files contain data from following DUs: {du_list}')
'''
# get some useful parameters from the data directory
data_dir_split = data_dir.split('/')
array = path_to_data_dir_split[4]
month = path_to_data_dir_split[5]
mode  = path_to_data_dir_split[6]
'''
######################################
# CORE FUNCTION TO SEARCH TRANSIENTS #
######################################

def std_search(trace, num_threshold, standard_separation): 
    # stop if there are no transients
    threshold = num_threshold * np.std(trace)
    exceed_threshold = np.abs(trace) > threshold
    if not np.any(exceed_threshold):
        return []
    
    # find trace positions where threshold is exceeded
    crossing_ids = np.flatnonzero(exceed_threshold)
    
    # find the separations between consecutive threshold crossings
    separations_crossings = np.diff(crossing_ids)

    # locate pulse indices in threshold crossing indices
    pulse_ids = np.flatnonzero(separations_crossings > standard_separation)
    pulse_ids = np.concatenate(([-1], pulse_ids, [len(crossing_ids)-1]))
    
    # preallocate the return list for time windows
    window_list = [[0, 0] for _ in range(len(pulse_ids) - 1)]

    # search all transients/pulses
    half_separation = standard_separation // 2
    for i in range(len(pulse_ids)-1):
        # get the start index of current pulse
        start_id = crossing_ids[pulse_ids[i]+1] - half_separation
        start_id = max(0, start_id) # fix the 1st pulse

        # get the stop index of current pulse
        stop_id = crossing_ids[pulse_ids[i+1]] + half_separation
        stop_id = min(len(trace)-1, stop_id) # fix the last pulse

        window_list[i] = [start_id, stop_id]

    return window_list

def noise_search(trace, noise, num_threshold, standard_separation): 
    # stop if there are no transients
    threshold = num_threshold * noise
    exceed_threshold = np.abs(trace) > threshold
    if not np.any(exceed_threshold):
        return []
    
    # find trace positions where threshold is exceeded
    crossing_ids = np.flatnonzero(exceed_threshold)
    
    # find the separations between consecutive threshold crossings
    separations_crossings = np.diff(crossing_ids)

    # locate pulse indices in threshold crossing indices
    pulse_ids = np.flatnonzero(separations_crossings > standard_separation)
    pulse_ids = np.concatenate(([-1], pulse_ids, [len(crossing_ids)-1]))
    
    # preallocate the return list for time windows
    window_list = [[0, 0] for _ in range(len(pulse_ids) - 1)]

    # search all transients/pulses
    half_separation = standard_separation // 2
    for i in range(len(pulse_ids)-1):
        # get the start index of current pulse
        start_id = crossing_ids[pulse_ids[i]+1] - half_separation
        start_id = max(0, start_id) # fix the 1st pulse

        # get the stop index of current pulse
        stop_id = crossing_ids[pulse_ids[i+1]] + half_separation
        stop_id = min(len(trace)-1, stop_id) # fix the last pulse

        window_list[i] = [start_id, stop_id]

    return window_list

##################################
# SEARCH TRANSIENTS IN ALL FILES #
##################################

# mask channels to disable channel 0 for GP13 data
mask_channel = np.array([False, True, True, True])

# make dictionaries to avoid looping through all files more than once
windows = {channel: {} for channel in channel_list}
gps_times = {}

# loop through all DUs to initiate dictionaries with empty lists
for du in du_list:
    for channel in channel_list:
        windows[channel][du] = []
    gps_times[du] = []

# loop through all files
for i, file in enumerate(file_list):
    # get info of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()

    # loop over all events in this file
    print(f'\n{i+1}/{num_files}: Looping over {num_entries} events/entries in {file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # assume that 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
    
        # get traces and search transients
        traces = np.array(tadc.trace_ch[0])[mask_channel]
        for i, channel in enumerate(channel_list):
            windows[channel][du].append(std_search(traces[i], num_threshold, standard_separation))
            #windows[channel][du].append(noise_search(traces[i], noises[i], num_threshold, standard_separation))
        gps_times[du].append(trawv.gps_time[0])

# loop through all DUs
for du in du_list:
    # save all info into a NPZ file
    npz_dir  = 'result2/'
    npz_file = 'DU{}_threshold{}_separation{}_date{}.npz'.format(du, num_threshold, standard_separation, specific_date[:8])
    np.savez(os.path.join(npz_dir, specific_date, npz_file),
             du = du,
             **{f'window{channel}': windows[channel][du] for channel in channel_list},
             gps_times = gps_times[du])
    print(f'\nSaved: {npz_dir}{specific_date}{npz_file}')

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")