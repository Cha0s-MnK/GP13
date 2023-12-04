###################
# IMPORT PACKAGES #
###################

from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime, time, timedelta
from glob import glob
import itertools
import json
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
np.set_printoptions(threshold=np.inf, precision=2, suppress=True, linewidth=128)
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
channels            = ['X', 'Y', 'Z'] # make dictionaries for easier indexing; X(north-south), Y(east-west), Z(up-down)
cutoff_frequency    = 50 # cut-off frequency of the high pass filter [MHz]
std_fluctuation     = 1.5 # tolerance of abnormal fluctuation for trace standard deviation [times]
linear_gain         = 10 # system's linear gain 
max_samples         = 512 # maximum number of samples in a transient/pulse
num_crossings       = 6 # least number of threshold crossings in a time window
num_run             = 94 # run number of GP13 data
num_samples         = 2048 # number of samples in a trace
num_threshold       = 4 # trigger threshold of the transient
sample_frequency    = 500 # sampling frequency in each trace [MHz]
standard_separation = 24 # required separation between closest threshold crossings in two pulses [#sample]
time_step           = 2 # time step of each sample [ns]

fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency # frequencies of the FFT [MHz]
time_axis     = np.arange(num_samples) * time_step # time axis of a trace

# mask channels to disable channel 0 for GP13 data
channel_mask = np.array([False, True, True, True])

# all/default dates
date_list_RUN92 = ['20231115', '20231116', '20231117', '20231118', '20231119', '20231120', '20231121', '20231122', '20231123', '20231124', '20231125', '20231126', '20231127', '20231128']
date_list_RUN93 = ['20231117', '20231118', '20231119', '20231120', '20231121', '20231122', '20231123', '20231124', '20231125', '20231126', '20231127', '20231128']
date_list_RUN94 = ['20231118', '20231121', '20231122', '20231123', '20231124', '20231125']

# all/default DUs
du_list_RUN92 = [1010, 1013, 1017, 1019, 1020, 1021, 1029, 1032, 1035, 1041, 1076, 1085]
du_list_RUN93 = [1076, 1085]
du_list_RUN94 = [1085]

# save directories
data_dir   = 'data'
result_dir = 'result'
plot_dir   = 'plot'

########################################
# DEFINE ARGUMENTS FOR SPECIFIC SCRIPT #
########################################

# get_trace.py & plot_trace.py
check_du_list   = [1076, 1085]
trace_wanted    = 'abnormal'
channel_mapping = {'X': 0, 'Y': 1, 'Z': 2}

# random_noise.py
num_sim = 16384

# get_noise.py
noise_data_dir   = f'data/RUN{num_run}'
noise_result_dir = f'result/noise'

# load the noises
with open(os.path.join(noise_result_dir, f'noise_RUN{num_run}.json'), 'r') as file:
    noises = json.load(file)

# search_transient.py
search_data_dir     = f'data/RUN{num_run}'
save_dir   = 'result/search/'
duration_result_dir = 'result/duration'

# plot_rate.py
rate_plot_files = 'result/search/*/*.npz'
rate_plot_dir   = 'plot/rate/'

# plot mean FFTs for each DU
fft_plot_dir = 'plot/fft/'

# get_fft.py
fft_result_dir  = 'result/fft'

# plot_fft.py
galaxy_dir   = 'data/galaxy'
galaxy_name  = ['VoutRMS2_NSgalaxy.npy', 'VoutRMS2_EWgalaxy.npy', 'VoutRMS2_Zgalaxy.npy']
fft_plot_dir = 'plot/fft'

###############################
# GET DUS FROM ROOT/NPZ FILES #
###############################

def get_root_dus(file):
    # enable GRANDLIB
    import grand.dataio.root_trees as rt

    # initiate TADC tree
    tadc = rt.TADC(file)

    # get the 1st entry
    tadc.get_entry(0)

    # get used DUs
    du_list = tadc.get_list_of_all_used_dus()
    return du_list

def get_roots_dus(file_list):
    # enable GRANDLIB
    import grand.dataio.root_trees as rt

    # create an empty set to store unique DUs
    du_set = set()

    # loop through all files to get used DUs
    for file in file_list:
        tadc = rt.TADC(file) # initiate TADC tree of this file
        tadc.get_entry(0) # get the entry from this file
        du_set.update(tadc.get_list_of_all_used_dus()) # update the set with all used DUs from this file

    # convert the set to a list in order
    du_list = sorted(list(du_set))

    # print and return used DUs
    print(f'\nROOT files contain data from following DUs: \n{du_list}\n')
    return du_list

def get_npz_dus(file_list):
    # create an empty set to store unique DUs
    du_set = set()

    # show all NPZ files and extract DUs
    print('\nPlot with data from following NPZ files:\n')
    for file in file_list:
        basename = os.path.basename(file)
        print(basename)
    
        # split the filename on underscore
        du = int(basename.split('_')[2][2:])
        du_set.add(du)

    # convert the set to a list in order
    du_list = sorted(list(du_set))

    # print the list of all used DUs
    print(f'\nNPZ files contain data from following DUs: \n{du_list}\n')

    # return DUs
    return du_list

###########################
# GET TIME AND CONVERT IT #
###########################

def get_root_DHtime(file):
    # assume that all ROOT filenames have the same pattern
    DHtime_str  = os.path.basename(file).split('test.')[1].split('.')[0]
    DHtime      = datetime.strptime(DHtime_str, '%Y%m%d%H%M%S')
    DHtime_flat = DHtime.strftime('%Y%m%d%H%M%S')
    return DHtime, DHtime_flat

def get_root_dates(file_list):
    # create an empty set to store unique dates
    date_set = set()

    for file in file_list:
        date_time, datetime_flat = get_root_datetime(file)
        date_set.add(datetime_flat[:8])

    # convert the set to a list in order
    date_list = sorted(list(date_set))

    # print and return the dates
    print(f'\nROOT files contain data from following dates: \n{date_list}\n')
    return date_list

def get_npz_DHtime(filename):
    # assume that all NPZ filenames have the same pattern
    basename    = os.path.basename(filename)
    DHtime_str  = basename.split('.npz')[0].split('_')[-1]
    DHtime      = datetime.strptime(DHtime_str, '%Y%m%d%H%M%S')
    DHtime_flat = DHtime.strftime('%Y%m%d%H%M%S')
    return DHtime, DHtime_flat

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
                     sample_frequency=sample_frequency, 
                     cutoff_frequency=cutoff_frequency):
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
                   std_fluctuation=std_fluctuation,
                   sample_frequency=sample_frequency, 
                   cutoff_frequency=cutoff_frequency):
    # apply the high-pass filter if filter is set to 'on'
    if filter == 'on':
        trace = high_pass_filter(trace=trace, 
                                 sample_frequency=sample_frequency, 
                                 cutoff_frequency=cutoff_frequency)

    # stop if there are no transients
    exceed_threshold = np.abs(trace) > threshold
    if np.sum(exceed_threshold) < num_crossings:
        return []

    # stop if this trace has an abnormal fluctuation of standard deviation
    if np.std(trace) > (std_fluctuation * threshold / num_threshold):
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

def search_windows_test(trace,
                        num_threshold=num_threshold,
                        filter_status='on',
                        standard_separation=standard_separation,
                        num_crossings=num_crossings,
                        std_fluctuation=std_fluctuation,
                        sample_frequency=sample_frequency, 
                        cutoff_frequency=cutoff_frequency):
    # apply the high-pass filter if filter is set to 'on'
    if filter_status == 'on':
        trace = high_pass_filter(trace=trace, 
                                 sample_frequency=sample_frequency, 
                                 cutoff_frequency=cutoff_frequency)

    # stop if there are no transients
    exceed_threshold = np.abs(trace) > num_threshold * np.std(trace)
    if np.sum(exceed_threshold) < num_crossings:
        return []

    # stop if this trace is abnormal/noisy
    if max(get_psd(trace)) > 1e-7:
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
        if np.sum(exceed_threshold[start_id:stop_id+1]) >= num_crossings:
            window_list.append([start_id, stop_id])

    return window_list

def rough_search_windows(trace,
                         threshold,
                         filter='on',
                         standard_separation=standard_separation,
                         num_crossings=num_crossings,
                         sample_frequency=sample_frequency, 
                         cutoff_frequency=cutoff_frequency):
    # apply the high-pass filter if filter is set to 'on'
    if filter == 'on':
        trace = high_pass_filter(trace=trace, 
                                 sample_frequency=sample_frequency, 
                                 cutoff_frequency=cutoff_frequency)

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

###############
# COMPUTE PSD #
###############

# get PSD of 1 trace
def get_psd(trace):
    # get FFTs
    fft = np.abs(np.fft.rfft(trace))

    # convert FFT from ADC units to Volts and correct for system linear gain
    fft = fft * adcu2v / linear_gain

    # compute power of FFT and normalize it to number of samples in a trace
    fft_power   = fft * fft / num_samples / num_samples

    # compute frequency bin width
    frequency_bin_width = fft_frequency[1] - fft_frequency[0]

    # normalize the power to frequency bin width to get the PSD
    fft_psd = fft_power / frequency_bin_width

    return fft_psd

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