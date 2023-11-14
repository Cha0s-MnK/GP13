###################
# IMPORT PACKAGES #
###################

from datetime import time
import glob
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
begin_time = wall_time.perf_counter()

####################
# DEFINE ARGUMENTS #
####################

num_crossings       = 2 # least number of threshold crossings in a time window
num_samples         = 2048 # number of samples in a trace
num_sim             = 10000  # number of simulations
num_threshold       = 3 # trigger threshold of the transient (times of noises)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)
time_axis           = np.arange(num_samples) * 2 # time axis of a trace
time_step           = 2 # time step of each sample (ns)

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
# SIMULATION OF A GAUSSIAN NOISE #
##################################

# generate a Gaussian noise
np.random.seed(0) # for reproducibility
noise_trace = np.random.normal(loc=0.0, scale=1.0, size=num_samples)

# calculate the threshold and search for time windows
threshold = num_threshold * np.std(noise_trace)
window_list = search_windows(noise_trace, threshold, standard_separation)

# create a figure with a specific size and plot the noise trace
fig, ax = plt.subplots(figsize=(36, 4), dpi=256)
ax.plot(time_axis, noise_trace, color='blue', label='noise trace')

# add horizontal dashed lines at +/-Threshold
ax.axhline(y=threshold, color='red', linestyle='--', label=f'+/-Threshold={threshold:.2f}')
ax.axhline(y=-threshold, color='red', linestyle='--')

# set X-limits
ax.set_xlim(-4, 4100)

# highlight the time windows in the time axis
for i, window in enumerate(window_list):
    start_id, stop_id = window # 1 sample = 2 nanoseconds
    start_time = 2 * start_id # 1 sample = 2 nanoseconds
    stop_time  = 2 * stop_id
    # label only the 1st axvspan
    if i == 0:
        ax.axvspan(start_time, stop_time, color='green', alpha=0.3, label='Time Window(s)')
    else:
        ax.axvspan(start_time, stop_time, color='green', alpha=0.3)

    # annotate time windows
    centre_time = (start_time + stop_time) / 2
    ylim = ax.get_ylim()
    ax.text(centre_time, 1.02*ylim[1], f'[{start_time:.1f}, {stop_time:.1f} ]', color='black', fontsize=10, ha='center', va='bottom')

# adjust the size of the x-tick and y-tick labels
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)

# set an x-tick every 250 units
ax.set_xticks(np.arange(min(time_axis), max(time_axis)+1, 250))

# add some comments
ax.legend(loc='upper right')
ax.set_xlabel('Time / ns', fontsize=16)
ax.set_ylabel('ADC Counts', fontsize=16)
ax.set_title(f'Gaussian Noise Trace with Detected Windows \n Threshold = {num_threshold}, Separation = {standard_separation}', fontsize=20, pad=20)

# adjust layout
plt.tight_layout()

# save the figure as a PNG file
plt.savefig('plot/random_noise.png')

######################################
# SIMULATION OF DIFFERENT THRESHOLDS #
######################################

# 1 simulation of a noise trace
def noise_sim(num_threshold):
    # generate a Gaussian noise trace
    noise_trace = np.random.normal(loc=0.0, scale=1.0, size=num_samples)

    # search time windows in the noise trace
    window_list = search_windows(noise_trace, num_threshold, standard_separation)

    return(len(window_list))

# main simulation function for all thresholds
def sim(threshold_list, num_sim, PNGname):
    transient_rate_list = []
    for num_threshold in threshold_list:
        average_windows = sum(noise_sim(num_threshold) for _ in range(num_sim)) / num_sim

        # convert average number of detected windows to transient rate
        transient_rate = average_windows / (num_samples * time_step) * 1e6 # GHz --> kHz

        print(f'\n When threshold={num_threshold} and separation={standard_separation},')
        print(f'\n the average transient rate over {num_sim} simulations is {transient_rate}.')
        transient_rate_list.append(transient_rate)

    # plot the relation between number of thresholds and transient rates
    plt.figure(dpi=256)
    plt.plot(threshold_list, transient_rate_list, marker='o', linestyle='-')
    plt.xlabel('Number of Threshold')
    plt.ylabel('Transient rate / kHz')
    plt.title(f'Relation between Threshold and Transient Rate \n Gaussian Noise, Separation = {standard_separation}')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'plot/{PNGname}.png')

# simulation 1
sim(np.arange(3.0, 7.5, 0.5), num_sim, 'relation1')

# simulation 2
num_sim *= 2
sim(np.arange(3.7, 6.6, 0.1), num_sim, 'relation2')
'''
# simulation 3
num_sim *= 3
sim(np.arange(4.2, 6.2, 0.05), num_sim, 'relation3')

# simulation 4
num_sim *= 5
sim(np.arange(4.74, 5.72, 0.02), num_sim, 'relation4')
'''
# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - begin_time
print(f"\nWhole program executed in: {run_time} seconds")