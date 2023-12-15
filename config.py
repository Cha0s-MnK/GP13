###################
# IMPORT PACKAGES #
###################

from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime, time, timedelta
from glob import glob
import itertools
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FixedLocator
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
num_samples         = 2048 # number of samples in a trace
sample_frequency    = 500 # sampling frequency in each trace [MHz]
time_step           = 2 # time step of each sample [ns]

fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency # frequencies of the FFT [MHz]
time_axis     = np.arange(num_samples) * time_step # time axis of a trace

######################
# RUNNING PARAMETERS #
######################

threshold1     = 81
threshold2     = 36
maxRMS         = 21
maxPSD         = 5e-9
num_cross      = 2+4 # minimum number of threshold crossings in a transient/pulse window
num_interval   = 16
num_separation = 750

run_list = [104] # run numbers of GP13 data
nums_run = '-'.join(str(num_run) for num_run in run_list)

wanted = 'pulse' # 'trivial' / 'pulse'

# save directories
data_dir   = 'data'
result_dir = 'result'
plot_dir   = 'plot'

# all/default dates
# date_list_RUN104 = ['20231205', '20231206', '20231207', '20231208', '20231209', '20231210']

# all/default DUs
# DU_list_RUN104 = [1010, 1013, 1017, 1019, 1020, 1021, 1029, 1031, 1032, 1035, 1041, 1075, 1085]
good_DU_list = [1010, 1013, 1017, 1019, 1020, 1021, 1029, 1032, 1035, 1041, 1085]
good_DU_list = [1010, 1013, 1017, 1019, 1020, 1021, 1032, 1041]

DU_styles = {
    1010: {'colour': 'green', 'linestyle': '-', 'offset': 120},
    1013: {'colour': 'red', 'linestyle': '--', 'offset': 110},
    1017: {'colour': 'pink', 'linestyle': '-.', 'offset': 100},
    1019: {'colour': 'black', 'linestyle': ':', 'offset': 90},
    1020: {'colour': 'orange', 'linestyle': '-', 'offset': 80},
    1021: {'colour': 'purple', 'linestyle': '--', 'offset': 70},
    1029: {'colour': 'cyan', 'linestyle': '-.', 'offset': 60},
    1031: {'colour': 'darkorange', 'linestyle': ':', 'offset': 50},
    1032: {'colour': 'brown', 'linestyle': '-', 'offset': 40},
    1035: {'colour': 'purple', 'linestyle': '--', 'offset': 30},
    1041: {'colour': 'gray', 'linestyle': '-.', 'offset': 20},
    1075: {'colour': 'yellow', 'linestyle': ':', 'offset': 10},
    1085: {'colour': 'blue', 'linestyle': '-', 'offset': 0}
}

galaxy_dir  = 'data/galaxy'
galaxy_name = ['VoutRMS2_NSgalaxy.npy', 'VoutRMS2_EWgalaxy.npy', 'VoutRMS2_Zgalaxy.npy']

# get_noise.py
#noise_data_dir   = f'data/RUN{num_run}'
#noise_result_dir = f'result/noise'

# load the noises
#with open(os.path.join(noise_result_dir, f'noise_RUN{num_run}.json'), 'r') as file:
#    noises = json.load(file)

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

def get_root_DUs(file):
    # enable GRANDLIB
    import grand.dataio.root_trees as rt

    # initiate TADC tree
    tadc = rt.TADC(file)

    # get the 1st entry
    tadc.get_entry(0)

    # get used DUs
    DU_list = tadc.get_list_of_all_used_DUs()
    return DU_list

def get_roots_DUs(file_list):
    import grand.dataio.root_trees as rt

    DU_set = set()  # use set to store unique DUs

    for file in file_list:
        tadc = rt.TADC(file) # initiate TADC tree of this file
        tadc.get_entry(0) # get the 1st entry from this file
        DU_set.update(tadc.get_list_of_all_used_DUs()) # update the set

    DU_list = sorted(list(DU_set))
    print(f'ROOT files contain data from following DUs: \n{DU_list}')

    return DU_list

def get_root_dates(file_list):
    date_set = set() # use set to store unique dates

    for file in file_list:
        basename = os.path.basename(file)
        CST_str   = basename.split('.root')[0].split('_')[1] # all ROOT filenames have the same pattern
        date_flat = CST_str[:8]
        date_set.add(date_flat)

    date_list = sorted(list(date_set))
    print(f'ROOT files contain data from following dates (UTC): \n{date_list}')

    return date_list

def get1entry_CST(trawv):
    GPStime = trawv.gps_time[0]
    CST     = GPStime2UTC(GPStime) + timedelta(hours=8)
    CSTflat = CST.strftime('%Y%m%d%H%M%S')
    return CST, CSTflat

#####################
# GET INFO FROM NPZ #
#####################

def get_npz_files(name, date_list=False, DU_list=False):
    file_list = []

    if date_list and DU_list:
        for num_run in run_list: # 'run_list' is always required
            for date in date_list:
                for DU in DU_list:
                    files = sorted(glob(os.path.join(result_dir, f'{name}', f'*RUN{num_run}_DU{DU}*{date}.npz')))
                    file_list.extend(files)
                    print(f'Load {len(files)} NPZ files of DU{DU} on {date} from RUN{num_run}.')

    elif date_list and DU_list == False:
        for num_run in run_list:
            for date in date_list:
                files = sorted(glob(os.path.join(result_dir, f'{name}', f'*RUN{num_run}*{date}.npz')))
                file_list.extend(files)
                print(f'Load {len(files)} NPZ files on {date} from RUN{num_run}.')
    
    elif date_list == False and DU_list:
        for num_run in run_list:
            for DU in DU_list:
                files = sorted(glob(os.path.join(result_dir, f'{name}', f'*RUN{num_run}_DU{DU}*.npz')))
                file_list.extend(files)
                print(f'Load {len(files)} NPZ files of DU{DU} from RUN{num_run}.')

    else:
        for num_run in run_list:
            files = sorted(glob(os.path.join(result_dir, f'{name}', f'*RUN{num_run}*.npz')))
            file_list.extend(files)
            print(f'Load {len(files)} NPZ files from RUN{num_run}.')

    return file_list

def get_npz_DU(file):
    basename = os.path.basename(file)
    DU = int(basename.split('_')[2][2:]) # all NPZ filenames have the same pattern

    return DU

def get_npz_DUs(file_list): 
    DU_set = set() # use set to store unique DUs

    for file in file_list:
        DU = get_npz_DU(file)
        DU_set.add(DU)

    DU_list = sorted(list(DU_set))
    print(f'NPZ files contain data from following DUs: \n{DU_list}')

    return DU_list

def get_npz_CST(file):
    basename = os.path.basename(file)
    CSTflat = basename.split('.npz')[0].split('_')[-1] # all NPZ filenames have the same pattern
    CST      = datetime.strptime(CSTflat, '%Y%m%d%H%M%S')
    return CST, CSTflat

def get_npz_dates(file_list): 
    date_set = set() # use set to store unique dates

    for file in file_list:
        basename  = os.path.basename(file)
        CSTflat  = basename.split('.npz')[0].split('_')[-1] # all NPZ filenames have the same pattern
        date_flat = CSTflat[:8]
        date_set.add(date_flat)
    
    date_list = sorted(list(date_set))
    print(f'NPZ files contain data from following dates: \n{date_list}')

    return date_list

###########################
# GET TIME AND CONVERT IT #
###########################

# convert GPS time to UTC
def GPStimes2UTCs(GPStimes):
    # vectorize the helper function
    return np.vectorize(GPStime2UTC)(GPStimes)

# a helper function
def GPStime2UTC(GPStime):
    # GPS time = UTC + 18s at present
    leap_seconds = 18
    return datetime.utcfromtimestamp(GPStime - leap_seconds)

###########################
# SEARCH FOR TIME WINDOWS #
###########################

def high_pass_filter(trace, sample_frequency=sample_frequency, cutoff_frequency=cutoff_frequency):
    # Nyquist-Shannon sampling theorem: maximum frequency that can be effectively sampled without aliasing when the signal is sampled at a given rate
    Nyquist_frequency = 0.5 * sample_frequency

    # design a Butterworth high-pass filter
    b, a = butter(4, cutoff_frequency / Nyquist_frequency, btype='high', analog=False)
    return filtfilt(b, a, trace)
'''
def search_window(trace,
                  filter_status='on',
                  threshold1=threshold1,
                  threshold2=threshold2,
                  num_crossing_min=num_crossing_min,
                  num_crossing_max=num_crossing_max,
                  num_interval=num_interval,
                  sample_frequency=sample_frequency,
                  cutoff_frequency=cutoff_frequency):
    if filter_status == 'on':
        trace = high_pass_filter(trace=trace, 
                                 sample_frequency=sample_frequency, 
                                 cutoff_frequency=cutoff_frequency)

    #
    if get1trace_RMS(trace) > 20 or max(get1trace_PSD(trace=trace)[205:615]) > 8e-9:
        return []

    # stop if there are no transients/pulses
    cross_threshold1 = np.abs(trace) > threshold1
    cross_threshold2 = np.abs(trace) > threshold2
    if np.sum(cross_threshold1) < 2 or np.sum(cross_threshold2) < num_crossing_min:
        return []
    
    # find trace positions where threshold is exceeded
    cross1_ids = np.flatnonzero(cross_threshold2)
    
    # find the separations between consecutive threshold crossings
    cross_separations = np.diff(cross1_ids)

    # locate pulse indices in threshold crossing indices
    pulse_ids = np.flatnonzero(cross_separations > num_interval)
    pulse_ids = np.concatenate(([-1], pulse_ids, [len(cross1_ids)-1]))

    window_list = []

    # search all transients/pulses
    for i in range(len(pulse_ids)-1):
        # get the start index of current pulse
        start_id = cross1_ids[pulse_ids[i]+1] - num_interval
        start_id = max(0, start_id) # fix the 1st pulse

        # get the stop index of current pulse
        stop_id = cross1_ids[pulse_ids[i+1]] + num_interval
        stop_id = min(len(trace)-1, stop_id) # fix the last pulse


        if np.sum(cross_threshold2[start_id:stop_id+1]) >= num_crossing_min:
            if np.sum(cross_threshold2[start_id:stop_id+1]) < num_crossing_max:
                if np.sum(cross_threshold1[start_id:stop_id+1]) >= 2:
                    #trigger_id = start_id + np.flatnonzero(cross_threshold1[start_id:stop_id+1])[0]
                    #if np.sum(cross_threshold2[start_id:trigger_id+1]) <= 2:
                    window_list.append([start_id, stop_id])
        
    return window_list
'''
def search_window_test(trace,
                  filter_status='on',
                  threshold1=threshold1,
                  threshold2=threshold2,
                  maxRMS=maxRMS,
                  num_cross=num_cross,
                  num_separation=num_separation,
                  num_interval=num_interval,
                  sample_frequency=sample_frequency,
                  cutoff_frequency=cutoff_frequency):
    if filter_status == 'on':
        trace = high_pass_filter(trace=trace, sample_frequency=sample_frequency, cutoff_frequency=cutoff_frequency)
 
    # stop if there are no transients/pulses
    cross_threshold1 = np.abs(trace) > threshold1
    cross_threshold2 = np.abs(trace) > threshold2
    if np.sum(cross_threshold1) < 2 or np.sum(cross_threshold2) < num_cross:
        return []

    # exclude abnormal traces
    if get1trace_RMS(trace=trace) > maxRMS or max(get1trace_PSD(trace=trace)) > maxPSD:
        return []
    
    # find trace positions where threshold1 is exceeded
    cross1_ids = np.flatnonzero(cross_threshold1)
    
    # find the separations between consecutive threshold crossings
    cross1_separations = np.diff(cross1_ids)

    # locate pulse indices in threshold crossing indices
    pulse_ids = np.flatnonzero(cross1_separations > num_separation)
    pulse_ids = np.concatenate(([-1], pulse_ids, [len(cross1_ids)-1]))

    window_list = []

    # search all transients/pulses
    for i in range(len(pulse_ids)-1):
        trigger_id = cross1_ids[pulse_ids[i]+1]
        if np.sum(cross_threshold2[trigger_id-num_interval:trigger_id]) > 2:
            continue

        # get the start index of current pulse
        start_id = trigger_id - num_interval
        start_id = max(0, start_id) # fix the 1st pulse

        # get the stop index of current pulse
        cross2_ids       = np.flatnonzero(cross_threshold2)[np.flatnonzero(cross_threshold2) >= trigger_id]
        cross2_intervals = np.diff(cross2_ids)
        pulse2_ids       = np.flatnonzero(cross2_intervals > num_interval)
        if len(pulse2_ids) == 0:
            stop_id = cross2_ids[-1] + num_interval
        else:
            stop_id = cross2_ids[pulse2_ids[0]] + num_interval
        stop_id = min(len(trace)-1, stop_id) # fix the last pulse

        if np.sum(cross_threshold1[start_id:stop_id+1]) >= 2 and np.sum(cross_threshold2[start_id:stop_id+1]) >= num_cross:
            window_list.append([start_id, stop_id])
        
    return window_list
#'''
###############
# COMPUTE PSD #
###############

# compute RMS of a data list
def get1trace_RMS(trace):
    return np.sqrt(np.mean(np.square(trace)))

# get PSD of 1 trace
def get1trace_PSD(trace):
    # get FFTs
    fft = np.abs(np.fft.rfft(trace))

    # convert FFT from ADC units to Volts and correct for system linear gain
    fft = fft * adcu2v / linear_gain

    # compute power of FFT and normalize it to number of samples in a trace
    fft_power   = fft * fft / num_samples / num_samples

    # compute frequency bin width
    frequency_bin_width = fft_frequency[1] - fft_frequency[0]

    # normalize the power to frequency bin width to get the PSD
    fft_PSD = fft_power / frequency_bin_width

    return fft_PSD

def load1day_RMS_PSD(file_list, DU_list, start_time, stop_time):
    # make dictionaries for easier indexing and initiate them
    RMSs_list        = {DU: {} for DU in DU_list}
    PSDs_list         = {DU: {} for DU in DU_list}
    times_list        = {DU: {} for DU in DU_list}
    temperatures_list = {DU: {} for DU in DU_list}
    mean_PSDs         = {DU: {} for DU in DU_list}
    for DU in DU_list:
        for channel in channels:
            RMSs_list[DU][channel]        = []
            PSDs_list[DU][channel]         = []
            times_list[DU][channel]        = []
            temperatures_list[DU][channel] = []
            mean_PSDs[DU][channel]         = []

    # get info from NPZ files
    for file in file_list:
        npz_file = np.load(file, allow_pickle=True)
        DU       = get_npz_DU(file)
        for channel in channels:
            RMSs_list[DU][channel]        = npz_file[f'RMSs_list{channel}']
            PSDs_list[DU][channel]         = npz_file[f'PSDs_list{channel}']
            times_list[DU][channel]        = npz_file[f'times_list{channel}']
            temperatures_list[DU][channel] = npz_file[f'temperatures_list{channel}']

            mean_PSDs[DU][channel]         = np.mean(np.array(PSDs_list[DU][channel]), axis=0)

    # use a mask to filter the data
    for file in file_list:
        for DU in DU_list:
            for channel in channels:
                time_mask                      = (times_list[DU][channel] >= start_time) & (times_list[DU][channel] < stop_time)
                times_list[DU][channel]        = times_list[DU][channel][time_mask]
                RMSs_list[DU][channel]        = RMSs_list[DU][channel][time_mask]
                PSDs_list[DU][channel]         = PSDs_list[DU][channel][time_mask]
                temperatures_list[DU][channel] = temperatures_list[DU][channel][time_mask]
    
    return times_list, RMSs_list, PSDs_list, temperatures_list, mean_PSDs

def get_filter_list(origin_list, PSDlist, RMSlist):
    RMSmask = np.array([RMS <= 20 for RMS in RMSlist])
    PSDmask = np.array([max(sublist[:410]) <= 8e-9 for sublist in PSDlist])
    mask    = RMSmask & PSDmask
    
    filter_list = np.array(origin_list)[mask].tolist()

    return filter_list

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
        print(f"Whole program has been executed in {run_time:.2f} seconds.")