# import necessary packages
import glob
import grand.dataio.root_trees as rt # GRANDLIB
import numpy as np
import os
import time as wall_time
start_time = wall_time.perf_counter()

#############################################
# SEARCH FOR TRANSIENTS IN GRAND DATA FILES #
#############################################

# GRANDLIB BRANCH: 'dev_io_root'
# This script reads out a GRAND data file in ROOT format. It searches for transient pulses in the time traces of the data 
# file, where each DU is treated separately. Here, a transient is defined as a pulse that exceeds 5 times the mean noise 
# level of the trace, called 5 sigma from this point onwards. For each trace, the number of 5 sigma transients is computed, as well as 
# the number of times this 5 sigma threshold is exceeded during the transient. The maximum peak position of the transient is
# also stored. All information is saved in an NPZ file for each DU.

####################
# DEFINE ARGUMENTS #
####################

# make dictionaries to avoid defining too many variables
channels  = ['X', 'Y', 'Z']

data_dir            = 'data/' # path of data directory containing ROOT files to analyze
specific_date       = '20231031/' # specify ROOT files of which date to analyze
num_threshold       = 5 # trigger threshold of the transient (times of noises)
noises              = [19.0, 19.0, 19.0] # noise level for 3 ADC channels (ADC counts)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)

# get ROOT files and their info to analyze
files     = sorted(glob.glob(os.path.join(data_dir, specific_date, '*.root')))
num_files = len(files)

# define the NPZ directory to store the computed PSDs
npz_dir = 'result2/'

# input DU list manually
du_list = [1010, 1013, 1017, 1019, 1020, 1021, 1029, 1031, 1032, 1035, 1041]
'''
#####################################
# GET ROOT FILES AND DU INFORMATION #
#####################################

# initiate TADC tree of the file list to get basic info
total_entries   = 0 # total number of entries across all files
max_dus         = 0 # maximum used DUs
id_max_dus_file = 0 # index of file that has maximum used DUs
for i, file in enumerate(files):
    tadc = rt.TADC(file)
    tadc.get_entry(0) # get the entry from current file
    total_entries += tadc.get_number_of_entries()
    cur_dus = len(tadc.get_list_of_all_used_dus()) # get used DUs from current file
    if cur_dus > max_dus:
        max_dus         = cur_dus
        id_max_dus_file = i

# get list of used DUs
tadc    = rt.TADC(files[id_max_dus_file])
du_list = tadc.get_list_of_all_used_dus()
print(f'\nROOT files contain data from following DUs: {du_list}')

# get some useful parameters from the data directory
path_to_data_dir_split = path_to_data_dir.split('/')
array = path_to_data_dir_split[4]
month = path_to_data_dir_split[5]
mode  = path_to_data_dir_split[6]
'''
######################################
# MAIN FUNCTION TO SEARCH TRANSIENTS #
######################################

def fast_search(trace, noise, num_threshold, standard_separation): 
    # stop if there are no transients
    threshold = num_threshold * noise
    exceed_threshold = np.abs(trace) > threshold
    if not np.any(exceed_threshold):
        return np.empty((0, 2), dtype=int)
    
    # find trace positions where threshold is exceeded
    crossing_ids = np.flatnonzero(exceed_threshold)
    
    # find the separations between consecutive threshold crossings
    separations_crossings = np.diff(crossing_ids)

    # locate pulse indices in threshold crossing indices
    pulse_ids             = np.flatnonzero(separations_crossings > standard_separation)
    pulse_ids             = np.concatenate(([-1], pulse_ids, [len(crossing_ids)-1]))
    
    # preallocate the return array for time windows
    windows = np.empty((len(pulse_ids)-1, 2), dtype=int)

    # search all transients/pulses
    half_separation = standard_separation // 2
    for i in range(len(pulse_ids)-1):
        # get the start index of current pulse
        start_id = crossing_ids[pulse_ids[i]+1] - half_separation
        start_id = max(0, start_id) # fix the 1st pulse

        # get the stop index of current pulse
        stop_id = crossing_ids[pulse_ids[i+1]] + half_separation
        stop_id = min(len(trace)-1, stop_id) # fix the last pulse

        windows[i] = [start_id, stop_id]

    return windows

#########################################
# COMPUTE AND PLOT MEAN PSD FOR EACH DU #
#########################################

# mask channels to disable channel 0 for GP13
mask_channel = np.array([False, True, True, True])

for du in du_list:
    # make dictionaries to avoid defining too many variables
    windows = {channel: [] for channel in channels}
    gps_times_list = []

    # loop over all files
    for file_id, file in enumerate(files):
        # get info of this file
        tadc        = rt.TADC(file)
        trawv       = rt.TRawVoltage(file)
        num_entries = tadc.get_number_of_entries()

        # loop over all events in file
        print(f'\n{file_id+1}/{num_files}: Looping over {num_entries} events in {file} for DU{du}')
        for entry in range(num_entries):
            # get info of this entry
            tadc.get_entry(entry)
            trawv.get_entry(entry)

            # assume that 1 entry corresponds to only 1 DU
            if du != tadc.du_id[0]:
                continue

            # get traces and search transients
            traces = np.array(tadc.trace_ch[0])[mask_channel]
            for i, channel in enumerate(channels):
                windows[channel].append(fast_search(traces[i], noises[i], num_threshold, standard_separation))
            gps_times_list.append(trawv.gps_time[0])
    
    # convert lists to numpy arrays
    for channel in channels:
        windows[channel] = np.array(windows[channel], dtype=object)
    gps_times = np.array(gps_times_list)

    # save total info in NPZ file
    npz_file = 'DU{}_threshold{}_separation{}_date{}.npz'.format(du, num_threshold, standard_separation, specific_date[:8])
    np.savez(os.path.join(npz_dir, npz_file),
             du = du,
             **{f'window{channel}': windows[channel] for channel in channels},
             #**{f'crossing{channel}': crossings[channel] for channel in channels},
             gps_times = gps_times)
    print(f'\nSaved NPZ file: {npz_dir}{npz_file}')

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")