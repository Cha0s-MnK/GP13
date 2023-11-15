###################
# IMPORT PACKAGES #
###################

from collections import defaultdict
from datetime import datetime, time, timedelta
import glob
import itertools
from matplotlib.lines import Line2D
# create custom legend items
custom_lines = [Line2D([0], [0], color='red', linestyle='dashed', lw=2),
                Line2D([0], [0], color='orange', linestyle='dashed', lw=2)]
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
plt.rc('font', size=10)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["axes.titlesize"] = 18
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.labelsize'] = 12
import numpy as np
np.set_printoptions(threshold=np.inf)
import os
from scipy.signal import butter, filtfilt
import time as wall_time
start_time = wall_time.perf_counter()

####################
# DEFINE ARGUMENTS #
####################

channels            = ['X', 'Y', 'Z'] # make dictionaries for easier indexing
cutoff_frequency    = 50 # cut-off frequency of the high pass filter; [MHz]
data_dir            = 'data/' # path of data directory containing ROOT files to analyze
npz_files           = 'result/*/*.npz' # path of all NPZ files to plot
num_crossings       = 3 # least number of threshold crossings in a time window
num_samples         = 2048 # number of samples in a trace
num_threshold       = 5 # trigger threshold of the transient
noises              = [25.0, 25.0, 55.0] # average noise level for 3 ADC channels; [ADC counts]
result_dir          = 'result/' # path of result directory containing NPZ files to plot
sample_frequency    = 500 # sampling frequency in each trace; [MHz]
standard_separation = 100 # required separation between closest threshold crossings in two pulses; [#sample]
time_step           = 2 # time step of each sample; [ns]

# mask channels to disable channel 0 for GP13 data
channel_mask = np.array([False, True, True, True])

# all/default dates
date_list = ['20231011/', '20231012/', '20231013/', '20231014/', '20231015/', '20231016/', '20231017/', '20231018/', 
             '20231019/', '20231020/', '20231027/', '20231028/', '20231029/', '20231030/', '20231031/']

# all/default DUs
du_list = ['1010', '1013', '1016', '1017', '1019', '1020', '1021', '1029', '1031', '1032', '1033', '1035', '1041']

####################################
# CORE FUNCTIONS TO SEARCH WINDOWS #
####################################

def high_pass_filter(trace):
    # Nyquist-Shannon sampling theorem: maximum frequency that can be effectively sampled without aliasing when the signal is sampled at a given rate
    Nyquist_frequency = 0.5 * sample_frequency

    # design a Butterworth high-pass filter
    b, a = butter(4, cutoff_frequency / Nyquist_frequency, btype='high', analog=False)
    return filtfilt(b, a, trace)

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

def iterative_search(trace, num_threshold, standard_separation): # not tested and used
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