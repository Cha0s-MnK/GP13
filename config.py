###################
# IMPORT PACKAGES #
###################

from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime, time, timedelta
from glob import glob
import itertools
import json
from matplotlib.colors import LogNorm
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

#####################
# DEFAULT ARGUMENTS #
#####################

adcu2v              = 1.8 / 16384 # convert from ADC units back to Volts
channels            = ['X', 'Y', 'Z'] # make dictionaries for easier indexing; X(north-south), Y(east-west), Z(up-down)
channel_mapping     = {'X': 0, 'Y': 1, 'Z': 2}
channel_mask        = np.array([False, True, True, True]) # mask channels to disable channel 0 for GP13 data
cutoff_frequency    = 50 # cut-off frequency of the high pass filter [MHz]
linear_gain         = 10 # linear gain of the system
num_crossings       = 6 # least number of threshold crossings in a time window
num_samples         = 2048 # number of samples in a trace
num_threshold       = 4 # trigger threshold of the transient
sample_frequency    = 500 # sampling frequency in each trace [MHz]
standard_separation = 24 # required separation between closest threshold crossings in two pulses [#sample]
time_step           = 2 # time step of each sample [ns]

fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency # frequencies of the FFT [MHz]
time_axis     = np.arange(num_samples) * time_step # time axis of a trace

######################
# RUNNING PARAMETERS #
######################

run_list = [104] # run numbers of GP13 data
nums_run = '-'.join(str(num_run) for num_run in run_list)

wanted = 'pulse'

du_list = [1076, 1085]

# save directories
data_dir   = 'data'
result_dir = 'result'
plot_dir   = 'plot'

# all/default dates
# date_list_RUN92 = ['20231117', '20231118', '20231119', '20231120', '20231121', '20231122', '20231123', '20231124', '20231125', '20231126', '20231127', '20231128']
# date_list_RUN93 = ['20231118', '20231119', '20231120', '20231121', '20231122', '20231123', '20231124', '20231125', '20231126', '20231127', '20231128']
# date_list_RUN94 = ['20231121', '20231122', '20231123', '20231124', '20231125']

# all/default DUs
# du_list_RUN92 = [1010, 1013, 1017, 1019, 1020, 1021, 1029, 1032, 1035, 1041]
# du_list_RUN93 = [1076]
# du_list_RUN94 = [1085]

########################################
# DEFINE ARGUMENTS FOR SPECIFIC SCRIPT #
########################################

# get_trace.py & plot_trace.py

# random_noise.py
num_sim = 16384

# get_noise.py
#noise_data_dir   = f'data/RUN{num_run}'
#noise_result_dir = f'result/noise'

# load the noises
#with open(os.path.join(noise_result_dir, f'noise_RUN{num_run}.json'), 'r') as file:
#    noises = json.load(file)

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

######################
# GET INFO FROM ROOT #
######################

def get_root_files(run_list=run_list, date_list=False):
    file_list = []

    if date_list:
        for num_run in run_list: # 'run_list' is always required
            for date in date_list:
                files = sorted(glob(os.path.join(data_dir, f'RUN{num_run}', f'GP13_{date}*RUN{num_run}*.root')))
                file_list.extend(files)
                print(f'Load {len(files)} ROOT files on {date} from RUN{num_run}.')

    else:
        for num_run in run_list:
            files = sorted(glob(os.path.join(data_dir, f'RUN{num_run}', f'GP13*RUN{num_run}*.root')))
            file_list.extend(files)
            print(f'Load {len(files)} ROOT files from RUN{num_run}.')

    return file_list

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
    import grand.dataio.root_trees as rt

    du_set = set()  # use set to store unique DUs

    for file in file_list:
        tadc = rt.TADC(file) # initiate TADC tree of this file
        tadc.get_entry(0) # get the 1st entry from this file
        du_set.update(tadc.get_list_of_all_used_dus()) # update the set

    du_list = sorted(list(du_set))
    print(f'ROOT files contain data from following DUs: \n{du_list}')

    return du_list

def get_root_dates(file_list):
    date_set = set() # use set to store unique dates

    for file in file_list:
        basename = os.path.basename(file)
        cst_str   = basename.split('.root')[0].split('_')[1] # all ROOT filenames have the same pattern
        date_flat = cst_str[:8]
        date_set.add(date_flat)

    date_list = sorted(list(date_set))
    print(f'ROOT files contain data from following dates (UTC): \n{date_list}')

    return date_list

# get GPS time for 1 entry and convert it to CST
def get1entry_cst(trawv):
    gps_time = trawv.gps_time[0]
    cst      = gps2utc1(gps_time) + timedelta(hours=8)
    cst_flat = cst.strftime('%Y%m%d%H%M%S')
    return cst, cst_flat

#####################
# GET INFO FROM NPZ #
#####################

def get_npz_files(name, date_list=False, du_list=False):
    file_list = []

    if date_list and du_list:
        for num_run in run_list: # 'run_list' is always required
            for date in date_list:
                for du in du_list:
                    files = sorted(glob(os.path.join(result_dir, f'{name}', f'{name}_RUN{num_run}_DU{du}*{date}.npz')))
                    file_list.extend(files)
                    print(f'Load {len(files)} NPZ files of DU{du} on {date} from RUN{num_run}.')

    elif date_list and du_list == False:
        for num_run in run_list:
            for date in date_list:
                files = sorted(glob(os.path.join(result_dir, f'{name}', f'{name}_RUN{num_run}*{date}.npz')))
                file_list.extend(files)
                print(f'Load {len(files)} NPZ files on {date} from RUN{num_run}.')
    
    elif date_list == False and du_list:
        for num_run in run_list:
            for du in du_list:
                files = sorted(glob(os.path.join(result_dir, f'{name}', f'{name}_RUN{num_run}_DU{du}*.npz')))
                file_list.extend(files)
                print(f'Load {len(files)} NPZ files of DU{du} from RUN{num_run}.')

    else:
        for num_run in run_list:
            files = sorted(glob(os.path.join(result_dir, f'{name}', f'{name}_RUN{num_run}*.npz')))
            file_list.extend(files)
            print(f'Load {len(files)} NPZ files from RUN{num_run}.')

    return file_list

def get_npz_du(file):
    basename = os.path.basename(file)
    du = int(basename.split('_')[2][2:]) # all NPZ filenames have the same pattern

    return du

def get_npz_dus(file_list): 
    du_set = set() # use set to store unique DUs

    for file in file_list:
        du = get_npz_du(file)
        du_set.add(du)

    du_list = sorted(list(du_set))
    print(f'NPZ files contain data from following DUs: \n{du_list}')

    return du_list

def get_npz_cst(file):
    basename = os.path.basename(file)
    cst_flat = basename.split('.npz')[0].split('_')[-1] # all NPZ filenames have the same pattern
    cst      = datetime.strptime(cst_flat, '%Y%m%d%H%M%S')
    return cst, cst_flat

def get_npz_dates(file_list): 
    date_set = set() # use set to store unique dates

    for file in file_list:
        basename  = os.path.basename(file)
        cst_flat  = basename.split('.npz')[0].split('_')[-1] # all NPZ filenames have the same pattern
        date_flat = cst_flat[:8]
        date_set.add(date_flat)
    
    date_list = sorted(list(date_set))
    print(f'NPZ files contain data from following dates: \n{date_list}')

    return date_list

###########################
# GET TIME AND CONVERT IT #
###########################

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

def search_windows_old(trace,
                   threshold,
                   filter='off',
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

def search_windows(trace,
                   num_threshold=num_threshold,
                   filter_status='on',
                   standard_separation=standard_separation,
                   num_crossings=num_crossings,
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

###############
# COMPUTE PSD #
###############

# compute RMS of a data list
def rms(data_list):
    return np.sqrt(np.mean(np.square(data_list)))

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