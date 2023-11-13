###################
# IMPORT PACKAGES #
###################

from collections import defaultdict
from datetime import datetime, time, timedelta
import glob
import grand.dataio.root_trees as rt # GRANDLIB
import itertools
import numpy as np
np.set_printoptions(threshold=np.inf)
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
plt.rc('font', size=10)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.labelsize"] = 14
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
import os
import time as wall_time
start_time = wall_time.perf_counter()

####################
# DEFINE ARGUMENTS #
####################

channel_list        = ['X', 'Y', 'Z'] # make dictionaries for easier indexing
data_dir            = 'data/20231027' # path of data directory containing ROOT files to check
specific_file       = 'GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231027122253.121_dat.root' # specify which ROOT file to check
num_samples         = 2048 # number of samples in a trace
num_threshold       = 5 # trigger threshold of the transient (times of noises)
noise_list          = [20.0, 20.0, 40.0] # average noise level for 3 ADC channels (ADC counts)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)ss
time_step           = 2 # time step of each sample (ns)

##############################
# GET ROOT FILES AND DU LIST #
##############################

# get ROOT files and their information to analyze
file_list = sorted(glob.glob(os.path.join(data_dir, specific_file)))
num_files = len(file_list)

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
print(f'\nROOT files contain data from following DUs: \n{du_list}\n')

#######################################
# CORE FUNCTIONS TO SEARCH TRANSIENTS #
#######################################

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

#####################################
# PLOT TRANSIENTS IN SELECTED FILES #
#####################################

# convert GPS times to UTCs; GPS time = UTC + 18s at present
def gps2utc(gps_times):
    gps2utc_v = np.vectorize(gps2utc1) # vectorize the helper function
    return gps2utc_v(gps_times)

# a helper function to convert single GPS time to UTC
def gps2utc1(gps_time):
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    return datetime.utcfromtimestamp(gps_time - leap_seconds)

# 
def plot_trace(traces, gps_time):
    # convert GPS time to UTC
    utc_time = gps2utc1(gps_time)

    # create a 3x1 layout for the subplots and adjust figure size
    fig, axes = plt.subplots(3, 1, figsize=(36, 12), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the data per channel
    for id, channel in enumerate(channel_list):
        # search time windows
        windows = std_search(traces[id], num_threshold, standard_separation)

        # locate corresponding axis
        ax = axes[id]

        # plot the trace
        ax.plot(time_axis, traces[id], color='blue', label='trace')

        # enable gird
        ax.grid(True)

        # add horizontal dashed lines at +/-Threshold of the trace
        threshold = num_threshold * np.std(traces[id])
        ax.axhline(y=threshold, color='red', linestyle='--', label=f'+/-Threshold={threshold:.2f}')
        ax.axhline(y=-threshold, color='red', linestyle='--')

        # highlight the time windows
        for id, window in enumerate(windows):
            start_id, stop_id = window
            # label only the 1st axvspan
            if id == 0:
                ax.axvspan(start_id, stop_id, color='green', alpha=0.3, label='Time Window(s)')
            else:
                ax.axvspan(start_id, stop_id, color='green', alpha=0.3)

            # annotate time windows
            centre_id = (start_id + stop_id) / 2
            ylim = ax.get_ylim()
            ax.text(centre_id, 1.02*ylim[1], f'[{start_id:.1f}, {stop_id:.1f} ]', color='black', fontsize=10, ha='center', va='bottom')
        
        # enable ticks on both left and right sides of the plot
        ax.tick_params(axis='y', labelleft=True, labelright=True)

        # set title on the right-hand side
        ax.text(1.01, 0.5, f'Channel {channel}', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=20, rotation=-90)
        
        # set X-limits
        ax.set_xlim(-4, 4100)

        # adjust the size of the x-tick and y-tick labels
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        # only set X-tick labels for the last subplot
        if id < len(channels) - 1:
            ax.set_xticklabels([])
        else:
            # set an x-tick every 100 units
            xticks = np.arange(min(time_axis), max(time_axis)+1, 250)
            ax.set_xticks(xticks)

        # add some comments
        ax.legend(loc='upper right')
        ax.set_xlabel('Time / ns', fontsize=16)
        ax.set_ylabel('ADC Counts', fontsize=16)
        ax.set_title(f'Gaussian Noise Trace with Detected Windows \n Threshold = {num_threshold}, Separation = {standard_separation}', fontsize=20, pad=20)

    # set X label and Y label for entire figure
    fig.text(0.5, 0.0, 'DunHuang Time in Date and Hour', ha='center', fontsize=20)
    fig.text(0.0, 0.5, 'Transient Rate / kHz', va='center', rotation='vertical', fontsize=20)
    
    # add the figure legend
    fig.legend(custom_lines, ['Sunrise', 'Sunset'], loc='upper right', fontsize=18, bbox_to_anchor=(0.98,1), bbox_transform=plt.gcf().transFigure)

    # add a main/super title for the entire figure
    plt.suptitle(f'Time Evolution of Transient Rates \n DU{du}, Threshold = {num_threshold}, Separation = {standard_separation}', fontsize=24)

    # adjust layout
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 1.0])

    # save the figure as a PNG file
    plt.savefig(f'plot/check/DU{du}_threshold{num_threshold}_separation{standard_separation}.png')
    print(f'Saved: plot/check/_DU{du}_threshold{num_threshold}_separation{standard_separation}.png')

    # close the figure to free up memory
    plt.close(fig)

# mask channels to disable channel 0 for GP13 data
mask_channel = np.array([False, True, True, True])

# make dictionaries to avoid looping through all files more than once
traces_list = {channel: {} for channel in channel_list}
windows_list = {channel: {} for channel in channel_list}
gps_times = {}

# loop through all DUs to initiate dictionaries with empty lists
for du in du_list:
    for channel in channel_list:
        traces_list[channel][du] = []
        windows_list[channel][du] = []
    gps_times[du] = []

# generate corresponding time axis
time_axis = np.arange(num_samples) * 2 # 1 sample = 2 nanoseconds

# loop through all files
for i, file in enumerate(file_list):
    # get info of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()

    # loop over all events in this file
    print(f'{i+1}/{num_files}: Looping over {num_entries} events/entries in {file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # assume that 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
    
        # get traces and search transients
        traces = np.array(tadc.trace_ch[0])[mask_channel]
        gps_time = trawv.gps_time[0]
        plot_trace(traces, gps_time)

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")