# system
import os
import glob
import argparse
import datetime

# scipy
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)
# grandlib
import grand.dataio.root_trees as rt

# time
import time
start_time = time.perf_counter()

##############################
# SET CUSTOM PLOT PARAMETERS #
##############################

plt.rcParams.update({'axes.labelsize': 30})
plt.rcParams.update({'xtick.labelsize': 30})
plt.rcParams.update({'ytick.labelsize': 30})
plt.rcParams.update({'axes.titlesize': 50})
plt.rcParams.update({'legend.fontsize': 30})

###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Search for transient in data traces and save information in NPZ files per DU.")

parser.add_argument('--data_dir',
                    dest='data_dir',
                    default='granddata/',
                    type=str,
                    help='Specify path of data directory containing ROOT files to analyze.')

parser.add_argument("--specific_files",
                    dest="specific_files",
                    default='GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231027122253.121_dat.root',
                    type=str,
                    help="Specify which ROOT files to analyze.")

parser.add_argument('--empty_channel',
                    dest='empty_channel',
                    default=0,
                    type=int,
                    help='Specify which of the 4 ADC channels is not branched.')

parser.add_argument('--num_threshold',
                    dest='num_threshold',
                    default=5,
                    type=int,
                    help='Specify the trigger threshold in units of sigma/STDs.')

parser.add_argument('--standard_separation',
                    dest='standard_separation',
                    default=100,
                    type=int,
                    help='Specify the required separation between two pulses in a trace in units of sample numbers.')       

parse_args = parser.parse_args()

data_dir            = parse_args.data_dir
specific_files      = parse_args.specific_files
empty_channel       = parse_args.empty_channel
num_threshold       = parse_args.num_threshold
standard_separation = parse_args.standard_separation

specific_date       = '202310271222'

#####################################
# GET ROOT FILES AND DU INFORMATION #
#####################################

# check if the data directory or the file list is empty
if not os.path.exists(data_dir):
    raise FileNotFoundError(f'\nData directory is not found: {data_dir}')
files = sorted(glob.glob(os.path.join(data_dir, specific_files)))

if not files:
    raise ValueError("\nProvided specific file list is empty. Please ensure that you provide a non-empty file list.")

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
tadc = rt.TADC(files[id_max_dus_file])
du_list = tadc.get_list_of_all_used_dus()
print(f'\nROOT files contain data from following DUs: {du_list}')

# create a channel mask to enable 3 ADC channels and disable channel 0 for GP13
mask_channel = np.ones(4, dtype=bool)
mask_channel[empty_channel] = False

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
    pulse_ids             = np.flatnonzero(separations_crossings > standard_separation)
    pulse_ids             = np.concatenate(([-1], pulse_ids, [len(crossing_ids)-1]))
    
    # preallocate the return list for time windows
    windows = [[0, 0] for _ in range(len(pulse_ids) - 1)]

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

# make dictionaries so that you don't have to loop over all files each time you want to analyze a different DU
channels  = ['X', 'Y', 'Z']
traces_all = {channel: {} for channel in channels}
windows = {channel: {} for channel in channels}
gps_times = {}
savedir = 'check_transients/10271222-5/'
num_files = len(files)

for file_id, file in enumerate(files):
    # get info of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()
    # loop over all events in this file
    print(f'\n{file_id+1}/{num_files}: Looping over {num_entries} events in {file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # assume that 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
        
        fig, axs = plt.subplots(len(channels), figsize=(40,20),  sharex=True)
        
        # get traces and search transients
        traces = np.array(tadc.trace_ch[0])[mask_channel]
        for i, channel in enumerate(channels):
            window = std_search(traces[i], num_threshold, standard_separation)
            trace = traces[i]
            threshold = num_threshold * np.std(trace)
            axs[i].plot(np.arange(2048),trace,label=f'{channel} trace')
            axs[i].axhline(y=threshold, color='r', linestyle='--', linewidth=2)
            axs[i].axhline(y=0-threshold, color='r', linestyle='--', linewidth=2)
            for w in window :
                if w != []:
                    start, end = w[0],w[1]
                    axs[i].axvspan(start, end, color='green', alpha=0.3)
            axs[i].set_ylabel('Amplitude')
            axs[i].legend()

        plt.xlabel('Sample Number')
        plt.suptitle(f'DU {du} - 12:22 - threshold {num_threshold} - {entry}',fontsize = 40)
        plt.savefig(f'{savedir}DU{du}_threshold{num_threshold}_{entry}.png')
        print(f'saved fig:{savedir}DU{du}_threshold{num_threshold}_{entry}.png')
        plt.close()
                    
        

