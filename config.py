###################
# IMPORT PACKAGES #
###################

from collections import defaultdict
from contextlib import contextmanager
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
from scipy.stats import norm
import time as wall_time
start_time = wall_time.perf_counter()

####################
# DEFINE ARGUMENTS #
####################

channels            = ['X', 'Y', 'Z'] # make dictionaries for easier indexing
cutoff_frequency    = 50 # cut-off frequency of the high pass filter; [MHz]
data_dir            = 'data/' # path of data directory containing ROOT files to analyze
filter_state        = 'off' # state of the high pass filter
max_samples         = 512 # maximum number of samples in a transient/pulse
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
du_list = [1010, 1013, 1016, 1017, 1019, 1020, 1021, 1029, 1031, 1032, 1033, 1035, 1041]

# good DUs
good_du_list = [1010, 1017, 1019, 1020, 1021, 1029, 1032, 1035]

# make a dictionary to store average noise level
noises = {channel: {} for channel in channels}

noises['X'][1010] = 20.93
noises['Y'][1010] = 29.84
noises['Z'][1010] = 53.79

noises['X'][1017] = 30.42
noises['Y'][1017] = 20.58
noises['Z'][1017] = 48.43

noises['X'][1019] = 31.25
noises['Y'][1019] = 24.55
noises['Z'][1019] = 67.11

noises['X'][1020] = 39.21
noises['Y'][1020] = 24.66
noises['Z'][1020] = 53.13

noises['X'][1021] = 33.70
noises['Y'][1021] = 32.39
noises['Z'][1021] = 88.26

noises['X'][1029] = 23.40
noises['Y'][1029] = 18.47
noises['Z'][1029] = 44.61

noises['X'][1032] = 42.08
noises['Y'][1032] = 40.57
noises['Z'][1032] = 105.35

noises['X'][1035] = 29.69
noises['Y'][1035] = 25.83
noises['Z'][1035] = 62.39

####################################
# DATA FILES AND SAVED DIRECTORIES #
####################################

# search for transients/pulses
data_dir          = 'data/'
search_result_dir = 'result/search/'

# plot transient/pulse rates
rate_npz_files  = 'result/rate/*/*.npz'
rate_plot_dir   = 'plot/rms/'

# get mean FFTs for each DU
#specific_file       = 'GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231028121653.161_dat.root'
specific_file       = 'GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231028001953.141_dat.root'

# analyze RMS
rms_data_files = 'data/20231012/*.root'
rms_result_dir = 'result/rms/20231012'
rms_npz_files  = 'result/rms/20231028/*.npz'
rms_plot_dir   = 'plot/rms/20231028/'

##############################
# GET ROOT FILES AND DU LIST #
##############################

def get_root_du(file_path):
    # enable GRANDLIB
    import grand.dataio.root_trees as rt
    
    # get ROOT files and their information
    file_list = sorted(glob.glob(file_path))
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

    # only consider good DUs
    du_list = [du for du in du_list if du in good_du_list]

    # print the list of all used good DUs
    print(f'\nROOT files contain data from following good DUs: \n{du_list}\n')

    # return ROOT files, number of ROOT files and DUs
    return file_list, num_files, du_list

####################################
# CORE FUNCTIONS TO SEARCH WINDOWS #
####################################

def high_pass_filter(trace, 
                     sample_frequency, 
                     cutoff_frequency):
    # Nyquist-Shannon sampling theorem: maximum frequency that can be effectively sampled without aliasing when the signal is sampled at a given rate
    Nyquist_frequency = 0.5 * sample_frequency

    # design a Butterworth high-pass filter
    b, a = butter(4, cutoff_frequency / Nyquist_frequency, btype='high', analog=False)
    return filtfilt(b, a, trace)

def search_windows(trace, 
                   threshold,
                   filter='off',
                   standard_separation=standard_separation,
                   num_crossings=num_crossings,
                   max_samples=max_samples,
                   sample_frequency=sample_frequency, 
                   cutoff_frequency=cutoff_frequency):
    # apply the high-pass filter if filter is set to 'on'
    if filter == 'on':
        trace = high_pass_filter(trace, sample_frequency, cutoff_frequency)

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

        # check if this time window contains enough crossings and does not exceed maximum number of samples
        if np.sum(exceed_threshold[start_id:stop_id+1]) >= num_crossings and (stop_id - start_id + 1) <= max_samples:
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

#######################
# RECORD RUNNING TIME #
#######################

@contextmanager
def record_run_time():
    start_time = wall_time.perf_counter()
    try:
        yield
    finally:
        end_time = wall_time.perf_counter()
        run_time = end_time - start_time
        print(f"\nWhole program has been executed in: {run_time:.2f} seconds")