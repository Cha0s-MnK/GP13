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
DURATION            = 4096 # constant; (ns)
num_threshold       = 5 # trigger threshold of the transient (times of noises)
standard_separation = 100 # required separation between two pulses in a trace (sample numbers)

# input DUs manually
du_list = ['1010', '1013', '1016', '1017', '1019', '1020', '1021', '1029', '1031', '1032', '1033', '1035', '1041']

# path of data directory containing NPZ files to plot
data_dir      = 'result2/*/'
specific_file = '*.npz'

# get all NPZ files
files = sorted(glob.glob(os.path.join(data_dir, specific_file)))

# show all NPZ files
print('\nData from following NPZ files:')
for file in files:
    print(os.path.basename(file))

################################################
# COMPUTE AND PLOT TRANSIENT RATES FOR EACH DU #
################################################

# convert GPS times to UTCs; GPS time = UTC + 18s at present
def gps2utc(gps_times):
    gps2utc_v = np.vectorize(gps2utc1) # vectorize the helper function
    return gps2utc_v(gps_times)

# a helper function to convert single GPS time to UTC
def gps2utc1(gps_time):
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    return datetime.utcfromtimestamp(gps_time - leap_seconds)

# compute DunHuang hours and transient rates
def compute_rate(utcs, windows_channel):
    # pair UTCs with time windows
    pairs_utc_window = [(utc, window) for utc, window in zip(utcs, windows_channel)]

    # make a dictionary to pair merge time windows with the same hour together
    dict_hour_window = defaultdict(list)
    for utc, window in pairs_utc_window:
        dict_hour_window[(utc.date(), utc.hour)].append(window)

    # build lists of UTC hours and transient rates
    utc_hours       = []
    transient_rates = []
    for (date, hour), windows in dict_hour_window.items():
        utc_hours.append(datetime.combine(date, time(hour=hour)))
        total_duration = len(windows) * DURATION
        windows        = list(itertools.chain(*windows)) # filter out empty time windows and reformat windows
        transient_rates.append(len(windows)/total_duration*1e6) # GHz --> kHz

    # convert UTCs to DunHuang local times
    DunHuang_hours = [hour + timedelta(hours=8) for hour in utc_hours]
    
    return DunHuang_hours, transient_rates

# plot transient rates for 3 ADC channels of 1 DU
def plot_rate(du, DunHuang_hours, transient_rates):
    # create a 3x1 layout for the subplots and adjust figure size
    fig, axes = plt.subplots(3, 1, figsize=(25, 10))

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the data per channel
    for id, channel in enumerate(channels):
        # reformat DunHuang hours and transient rates to plot
        DunHuang_hours[channel][du] = list(itertools.chain(*DunHuang_hours[channel][du]))
        transient_rates[channel][du] = list(itertools.chain(*transient_rates[channel][du]))

        # locate corresponding axis and plot
        ax = axes[id]
        ax.scatter(DunHuang_hours[channel][du], transient_rates[channel][du])
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H')) # set the date format
        ax.set_xlabel('DunHuang Time in date and hour')
        ax.set_ylabel('Transient Rate / kHz')
        ax.set_title(f'Channel {channel}')
        ax.grid(True)

        # rotate X labels
        for tick in ax.get_xticklabels():
            tick.set_rotation(15)
        '''
        # add vertical lines at sunrise and sunset
        ylim = ax.get_ylim()
        ax.vlines([datetime(2023, 10, 27, 8, 00),
                   datetime(2023, 10, 27, 19, 00),
                   datetime(2023, 10, 28, 8, 00),
                   datetime(2023, 10, 28, 19, 00)], ymin=ylim[0], ymax=ylim[1], color = 'r', linestyles='dashed')
        '''
    # adjust the spacing between the subplots to accommodate the rotated X labels
    plt.subplots_adjust(hspace=0.4)

    # add a main/super title for the entire figure
    plt.suptitle(f'Transient Rate Evolution of DU{du} with {num_threshold} Threshold {standard_separation} Separation', fontsize=18)

    # save the figure as a PNG file
    plt.savefig(f'plot2/rate_DU{du}_threshold{num_threshold}_separation{standard_separation}.png')

    # close the figure to free up memory
    plt.close(fig)

# make dictionaries for easier indexing
channels = ['X', 'Y', 'Z']
DunHuang_hours = {channel: {} for channel in channels}
transient_rates = {channel: {} for channel in channels}
gps_times = {}

# loop through all DUs to initiate dictionaries with empty lists
for du in du_list:
    for channel in channels:
        DunHuang_hours[channel][du] = []
        transient_rates[channel][du] = []
    gps_times[du] = []

# loop through all files to compute transient rates
for file in files:
    # get information from this file
    npz_file = np.load(file, allow_pickle=True)
    du = npz_file['du']
    windows = {channel: npz_file[f'window{channel}'] for channel in channels}
    gps_times = npz_file['gps_times']

    # convert GPS times to UTCs
    utcs = gps2utc(gps_times)

    # compute transient rates for 3 ADC channels
    for channel in channels:
        hours, rates = compute_rate(utcs, windows[channel])
        DunHuang_hours[channel][du].append(hours)
        transient_rates[channel][du].append(rates)
    
# loop through all DUs to plot transient rates for 3 ADC channels
for du in du_list:
    plot_rate(du, DunHuang_hours, transient_rates)

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")