# packages
from collections import defaultdict
from datetime import datetime, time, timedelta
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
start_time = wall_time.perf_counter()

# constants and parameters
num_threshold       = 3 # trigger threshold of the transient (times of noises)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)

# generate a Gaussian noise
np.random.seed(0) # for reproducibility
num_samples = 2048
noise_trace = np.random.normal(loc=0.0, scale=1.0, size=num_samples)

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

# search noise windows in the noise trace
window_list = std_search(noise_trace, num_threshold, standard_separation)

# generate corresponding time axis
time_axis = np.arange(num_samples) * 2 # 1 sample = 2 nanoseconds

# create a figure with a specific size and plot the noise trace
fig, ax = plt.subplots(figsize=(32, 4))
ax.plot(time_axis, noise_trace, color='blue')

# calculate the standard deviation of the noise trace and add horizontal dashed lines at +/-STD
std = np.std(noise_trace)
ax.axhline(y=std, color='red', linestyle='--', label=f'+/-STD={std:.2f}')
ax.axhline(y=-std, color='red', linestyle='--')
ax.legend(loc='upper right')

# set X-limits
ax.set_xlim(-4, 4100)

# highlight the time windows
for window in window_list:
    start_id, stop_id = window
    ax.axvspan(start_id, stop_id, color='green', alpha=0.3)

    # annotate time windows
    centre_id = (start_id + stop_id) / 2
    ylim = ax.get_ylim()
    ax.text(centre_id, 1.02*ylim[1], f'[{start_id:.1f}, {stop_id:.1f} ]', color='black', fontsize=10, ha='center', va='bottom')

# adjust the size of the x-tick and y-tick labels
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)

# set an x-tick every 100 units
xticks = np.arange(min(time_axis), max(time_axis)+1, 250)
ax.set_xticks(xticks)

# add comments
ax.set_xlabel('Time / ns', fontsize=16)
ax.set_ylabel('ADC Counts', fontsize=16)
ax.set_title(f'Gaussian Noise Trace with Detected Windows \n Threshold = {num_threshold}, Separation = {standard_separation}', fontsize=20, pad=20)

# adjust layout
plt.tight_layout()

# save the figure as a PNG file
plt.savefig('plot2/random_noise.png')

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")