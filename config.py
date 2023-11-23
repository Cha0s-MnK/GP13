###################
# IMPORT PACKAGES #
###################

from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime, time, timedelta
from glob import glob
import itertools
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["axes.titlesize"] = 18
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.labelsize'] = 14
import numpy as np
np.set_printoptions(threshold=np.inf)
import os
import re
from scipy.signal import butter, filtfilt
from scipy.stats import norm
import time as wall_time
start_time = wall_time.perf_counter()

####################
# DEFINE ARGUMENTS #
####################

adcu2v              = 1.8 / 16384 # convert from ADC units back to Volts
channels            = ['X', 'Y', 'Z'] # make dictionaries for easier indexing
cutoff_frequency    = 50 # cut-off frequency of the high pass filter [MHz]
fluctuation         = 1.5 # tolerance of abnormal fluctuation
linear_gain         = 10 # system's linear gain 
max_samples         = 512 # maximum number of samples in a transient/pulse
num_crossings       = 2 # least number of threshold crossings in a time window
num_samples         = 2048 # number of samples in a trace
num_threshold       = 5 # trigger threshold of the transient
sample_frequency    = 500 # sampling frequency in each trace [MHz]
standard_separation = 50 # required separation between closest threshold crossings in two pulses [#sample]
time_step           = 2 # time step of each sample [ns]

fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency # frequencies of the FFT [MHz]
time_axis     = np.arange(num_samples) * time_step # time axis of a trace

# mask channels to disable channel 0 for GP13 data
channel_mask = np.array([False, True, True, True])

# all/default dates
date_list = ['20231011', '20231012', '20231013', '20231014', '20231015', '20231016', '20231017', '20231018', 
             '20231019', '20231020', '20231021', '20231022', '20231027', '20231028', '20231029', '20231030', 
             '20231031', '20231115', '20231116', '20231117', '20231119', '20231120', '20231121']

# all/default DUs
du_list = [1010, 1013, 1016, 1017, 1019, 1020, 1021, 1029, 1031, 1032, 1033, 1035, 1041, 1076, 1085]

# good DUs
good_du_list = [1010, 1013, 1016, 1017, 1019, 1020, 1021, 1029, 1031, 1032, 1033, 1035, 1041, 1076, 1085]

# make a dictionary to store average noise level
noisesX = {du: 30 for du in good_du_list}
noisesY = {du: 5 for du in good_du_list}
noisesZ = {du: 5 for du in good_du_list}
noises = {'X': noisesX, 'Y': noisesY, 'Z': noisesZ}

########################################
# DEFINE ARGUMENTS FOR SPECIFIC SCRIPT #
########################################

# search.py
search_data_dir   = 'data/'
search_result_dir = 'result/search/'

# plot_rate.py
custom_sun      = [Line2D([0], [0], color='red', linestyle='dashed', lw=2), # create custom legend items
                   Line2D([0], [0], color='orange', linestyle='dashed', lw=2)]
rate_plot_files = 'result/search/*/*.npz'
rate_plot_dir   = 'plot/rate/'

# plot mean FFTs for each DU
fft_plot_dir = 'plot/fft/'

# check.py
good_du          = 1010
check_data_files = 'data/20231120/*.root'
check_result_dir = 'result/check/20231120/'

# plot_check.py
check_plot_dir   = 'plot/check/20231120/'

# get_fft.py
fft_result_dir  = 'result/fft'

# plot_fft2.py
galaxy_sim_dir  = 'data/galaxy'
galaxy_sim_name = ['VoutRMS2_NSgalaxy.npy', 'VoutRMS2_EWgalaxy.npy', 'VoutRMS2_Zgalaxy.npy']
fft_plot_dir    = 'plot/fft'

# get_amplitude.py
amplitude_result_dir = 'result/amplitude'

# analyze RMS
rms_data_files = 'data/20231012/*.root'
rms_result_dir = 'result/rms/20231012'
rms_npz_files  = 'result/rms/20231028/*.npz'
rms_plot_dir   = 'plot/rms/20231028/'

#################################
# GET DUS FROM ROOT/NPZ FILES #
#################################

def get_root_du(files):
    # enable GRANDLIB
    import grand.dataio.root_trees as rt

    # get ROOT files
    file_list = sorted(glob(files))

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

    # return DUs
    return du_list

def get_npz_du(files):
    # get NPZ files
    file_list = sorted(glob(files))

    # create an empty set to store unique DUs
    du_set = set()

    # show all NPZ files and extract DUs
    print('\nPlot with data from following NPZ files:\n')
    for file in file_list:
        basename = os.path.basename(file)
        print(basename)
    
        # split the filename on underscore
        du = int(basename.split('_')[1][2:])
        du_set.add(du)

    # convert the set to a list in order
    du_list = sorted(list(du_set))

    # only consider good DUs
    du_list = [du for du in du_list if du in good_du_list]

    # print the list of all used DUs
    print(f'\nNPZ files contain data from following DUs: \n{du_list}\n')

    # return DUs
    return du_list

###########################
# GET TIME AND CONVERT IT #
###########################

def get_root_datetime(filename):
    # assume all ROOT filenames have the same pattern
    datetime_str  = filename.split('test.')[1].split('.')[0]
    date_time     = datetime.strptime(datetime_str, '%Y%m%d%H%M%S')
    datetime_flat = date_time.strftime('%Y%m%d%H%M%S')
    return date_time, datetime_flat

def get_npz_datetime(basename):
    # assume all NPZ filenames have the same pattern
    datetime_str  = basename.split('.')[0].split('_')[-1]
    date_time     = datetime.strptime(datetime_str, '%Y%m%d%H%M%S')
    datetime_flat = date_time.strftime('%Y%m%d%H%M%S')
    return date_time, datetime_flat

# convert GPS time to UTC
def gps2utc(gps_times):
    # vectorize the helper function
    gps2utc_v = np.vectorize(gps2utc1)
    return gps2utc_v(gps_times)

# a helper function
def gps2utc1(gps_time):
    # GPS time = UTC + 18s at present
    leap_seconds = 18
    return datetime.utcfromtimestamp(gps_time - leap_seconds)

###########################
# SEARCH FOR TIME WINDOWS #
###########################

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

    # stop if this is a bad trace
    if np.std(trace) > fluctuation * threshold / num_threshold:
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