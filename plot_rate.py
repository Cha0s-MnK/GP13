# import necessary packages
import argparse
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

DURATION = 4096 # constant; (ns)

###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Load NPZ files and plot corresponding transient rates in 3 ADC channels.")

parser.add_argument('--data_dir',
                    dest='data_dir',
                    default='result/',
                    type=str,
                    help='Specify the path of the directory containing the NPZ files to analyze.')

parser.add_argument("--specific_files",
                    dest="specific_files",
                    default='*.npz',
                    type=str,
                    help="Specify which NPZ files to analyze.")

parse_args = parser.parse_args()

data_dir       = parse_args.data_dir
specific_files = parse_args.specific_files

#####################################
# GET ROOT FILES AND DU INFORMATION #
#####################################

# check if the data directory or the file list is empty
if not os.path.exists(data_dir):
    raise FileNotFoundError(f'\nData directory is not found: {data_dir}')
files = sorted(glob.glob(os.path.join(data_dir, specific_files)))
if not files:
    raise ValueError("\nProvided specific file list is empty. Please ensure that you provide a non-empty file list.")

################################################
# COMPUTE AND PLOT TRANSIENT RATES FOR EACH DU #
################################################

# parse DU, num_threshold, and standard_separation information from filename
def get_file_info(filename):
    filename_parts = os.path.basename(filename).split('_')
    du = filename_parts[0][2:]
    threshold = filename_parts[1][9:]
    separation = filename_parts[2][10:].split('.')[0]
    return du, threshold, separation

def gps2utc(gps_times): # convert GPS times to UTCs; GPS time = UTC + 18s at present
    gps2utc_v = np.vectorize(gps2utc1) # vectorize the helper function
    return gps2utc_v(gps_times)

def gps2utc1(gps_time): # a helper function to convert single GPS time to UTC
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    return datetime.utcfromtimestamp(gps_time - leap_seconds)

def plot_rate(du, threshold, separation, utcs, sup_windows): # compute and plot transient rates for 3 ADC channels of selected DU
    # create a 2x2 layout for the subplots and adjust figure size as needed
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # remove the last/4th ax as we only have 3 ADC channels to plot
    fig.delaxes(axes[1,1])

    # flatten axes array for easier indexing
    axes = axes.flatten()

    for id, channel in enumerate(channels):
        ax = axes[id]

        # use a list comprehension to filter out empty lists and pair time with non-empty time windows
        windows_channel = sup_windows[channel]
        pairs_utc_window = [(utc, window) for utc, window in zip(utcs, windows_channel)]

        # make a dictionary to pair hours with corresponding time windows and durations
        dict_hour_window = defaultdict(list)
        for utc, window in pairs_utc_window:
            dict_hour_window[(utc.date(), utc.hour)].append(window)

        # create lists of hours and transient rates
        date_hour      = []
        transient_rate = []
        for (date, hour), windows in dict_hour_window.items():
            date_hour.append(datetime.combine(date, time(hour=hour)))
            total_duration = len(windows) * DURATION
            windows        = list(itertools.chain(*windows)) # merge all time windows in the same time unit
            transient_rate.append(len(windows)/total_duration*1e6) # (10^9 Hz) --> kHz

        # plot the data per channel
        DunHuang_hour = [hour + timedelta(hours=8) for hour in date_hour]
        ax.plot(DunHuang_hour, transient_rate)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H')) # set the date format
        ax.set_xlabel('Time by date and hour')
        ax.set_ylabel('Transient Rate / kHz')
        ax.set_title(f'Channel {channel}')
        ax.grid(True)

        # rotate X labels
        for tick in ax.get_xticklabels():
            tick.set_rotation(15)

        # add vertical lines at sunrise and sunset
        ylim = ax.get_ylim()
        ax.vlines([datetime(2023, 10, 26, 19, 00),
                   datetime(2023, 10, 27, 8, 00),
                   datetime(2023, 10, 27, 19, 00),
                   datetime(2023, 10, 28, 8, 00),
                   datetime(2023, 10, 28, 19, 00)], ymin=ylim[0], ymax=ylim[1], color = 'r', linestyles='dashed')

    # adjust the spacing between the subplots to accommodate the rotated X labels
    plt.subplots_adjust(hspace=0.4)

    # add a main/super title for the entire figure
    plt.suptitle(f'Transient Rate Evolution of DU{du} with {threshold} Threshold {separation} Separation', fontsize=18)

    # save the figure as a PNG file
    plt.savefig(f'result/transient_rate_DU{du}_threshold{threshold}_separation{separation}.png')

    # close the figure to free up memory
    plt.close(fig)

# make dictionaries
channels = ['X', 'Y', 'Z']

# loop through files
for file in files:
    # extract information
    du, threshold, separation = get_file_info(file)
    npz_file = np.load(file, allow_pickle=True)
    windows = {channel: npz_file[f'window{channel}'] for channel in channels}
    gps_times = npz_file['gps_times']

    # convert GPS times to UTCs
    utcs = gps2utc(gps_times)

    # plot transient rates for 3 ADC channels
    plot_rate(du, threshold, separation, utcs, windows)

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")

'''
# use a list comprehension to filter out empty lists
window_x = np.array([a_list for a_list in npz_file['window_x'] if len(a_list) > 0])
window_y = np.array([a_list for a_list in npz_file['window_y'] if len(a_list) > 0])
window_z = np.array([a_list for a_list in npz_file['window_z'] if len(a_list) > 0])
trace_window_x = np.array([a_list for a_list in npz_file['trace_window_x'] if len(a_list) > 0])
trace_window_y = np.array([a_list for a_list in npz_file['trace_window_y'] if len(a_list) > 0])
trace_window_z = np.array([a_list for a_list in npz_file['trace_window_z'] if len(a_list) > 0])
start_x, end_x = int(window_x[0][0][0]), int(window_x[0][0][1])
start_y, end_y = int(window_y[0][0][0]), int(window_y[0][0][1])
start_z, end_z = int(window_z[0][0][0]), int(window_z[0][0][1])
list_window_x = list(range(start_x, end_x + 1))
list_window_y = list(range(start_y, end_y + 1))
list_window_z = list(range(start_z, end_z + 1))

# create a figure and a set of subplots
fig, ax = plt.subplots()

# Plot y versus x as lines and/or markers
ax.plot(list_window_x, trace_window_x[0][0])

# set the title and axis labels
ax.set_title('Transients')
ax.set_xlabel('Sample number')
ax.set_ylabel('Signal')

# save the figure
plt.savefig('result/test.png')
'''