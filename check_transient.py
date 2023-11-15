# system
import os
import glob
import argparse
from datetime import datetime, time, timedelta

# scipy
from scipy.signal import butter, filtfilt
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
plt.rcParams.update({'axes.titlesize': 30})
plt.rcParams.update({'legend.fontsize': 30})

####################
# DEFINE ARGUMENTS #
####################

channels            = ['X', 'Y', 'Z'] # make dictionaries for easier indexing
data_dir            = 'granddata/' # path of data directory containing ROOT files to analyze
specific_files      = 'GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231027122253.121_dat.root' #specific file
specific_date       = specific_files.split('.')[-3][4:12]
empty_channel       = 0
num_crossings       = 2 # least number of threshold crossings in a time window
num_threshold       = 5 # trigger threshold of the transient (times of noises)
noises              = [20.0, 20.0, 40.0] # average noise level for 3 ADC channels (ADC counts)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)

cutoff_frequency = 30  #filter out what's below, in MHz
sampling_rate = 500 #sample every 2 ns


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

#########################################
# COMPUTE AND PLOT MEAN PSD FOR EACH DU #
#########################################
# convert GPS times to UTCs; GPS time = UTC + 18s at present
def gps2utc(gps_times):
    gps2utc_v = np.vectorize(gps2utc1) # vectorize the helper function
    return gps2utc_v(gps_times)

# a helper function to convert single GPS time to UTC
def gps2utc1(gps_time):
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    return datetime.utcfromtimestamp(gps_time - leap_seconds)
    
def high_pass_filter(trace, cutoff_frequency, sampling_rate):
    nyquist = 0.5 * sampling_rate
    normal_cutoff = cutoff_frequency / nyquist
    b, a = butter(4, normal_cutoff, btype='high', analog=False)
    filtered_trace = filtfilt(b, a, trace)
    return filtered_trace


# make dictionaries so that you don't have to loop over all files each time you want to analyze a different DU
channels  = ['X', 'Y', 'Z']
traces_all = {channel: {} for channel in channels}
windows = {channel: {} for channel in channels}
gps_times = {}
savedir = f'check_transients/{specific_date}-{num_threshold}-filter{cutoff_frequency}/'
num_files = len(files)

for file_id, file in enumerate(files):
    # get info of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()
    # loop over all events in this file
    print(f'\n{file_id+1}/{num_files}: Looping over {num_entries} events in {file}')
    for entry in range(100):
        # get info of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # assume that 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
        
        fig, axs = plt.subplots(len(channels), figsize=(40,20),  sharex=True)
        
        #get time
        gps_times = trawv.gps_time[0]
        utctime = gps2utc(gps_times)
        hourdiff = timedelta(hours=8)
        localtime = utctime + hourdiff
        
        
        # get traces and search transients
        traces = np.array(tadc.trace_ch[0])[mask_channel]
        for i, channel in enumerate(channels):
            filtered_trace = high_pass_filter(traces[i], cutoff_frequency, sampling_rate)
            window = search_windows(filtered_trace, num_threshold*np.std(filtered_trace), standard_separation)
            threshold = num_threshold * np.std(filtered_trace)
            
            axs[i].plot(np.arange(2048),filtered_trace)
            axs[i].axhline(y=threshold, color='r', linestyle='--', linewidth=2, label=f'Threshold = {threshold:.2f}')
            axs[i].axhline(y=0-threshold, color='r', linestyle='--', linewidth=2)
            for w in window :
                if w != []:
                    start, end = w[0],w[1]
                    axs[i].axvspan(start, end, color='green', alpha=0.3)
            axs[i].set_ylabel('Amplitude')
            axs[i].legend(loc='upper right')
            axs[i].set_title(f'{channel} channel')

        plt.xlabel('Sample Number')
        plt.suptitle(f'DU {du} - {localtime.hour}:{localtime.minute}:{localtime.second} - threshold {num_threshold} - filter{cutoff_frequency}',fontsize = 40)
        plt.savefig(f'{savedir}DU{du}_threshold{num_threshold}_{localtime.hour}:{localtime.minute}:{localtime.second} - filter{cutoff_frequency}.png')
        print(f'saved fig:{savedir}DU{du}_threshold{num_threshold}_{localtime.hour}:{localtime.minute}:{localtime.second}_filter{cutoff_frequency}.png')
        plt.close()




