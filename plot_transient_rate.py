# import necessary packages
import argparse
from collections import defaultdict
from datetime import datetime, time
import glob
import grand.dataio.root_trees as rt # GRANDLIB
import itertools
import numpy as np
np.set_printoptions(threshold=np.inf)
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
plt.rc('font', size=10)
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
import os
import time as wall_time
start_time = wall_time.perf_counter()

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

# get list of used DUs
du_list = [file[9:13] for file in files]
print(f'\nNPZ files contain data from following DUs: {du_list}')

################################################
# COMPUTE AND PLOT TRANSIENT RATES FOR EACH DU #
################################################

def gps2utc(gps_times): # convert GPS times to UTCs; GPS time = UTC + 18s at present
    gps2utc_v = np.vectorize(gps2utc1) # vectorize the helper function
    return gps2utc_v(gps_times)

def gps2utc1(gps_time): # a helper function to convert single GPS time to UTC
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    return datetime.utcfromtimestamp(gps_time - leap_seconds)

def plot_transient_rate(utcs, windows, channel): # compute and plot transient rate for 3 ADC channels of selected DU
    # use a list comprehension to filter out empty lists and pair time with non-empty time windows
    windows_channel = windows[channel]
    pairs_utc_window = [(utc, window) for utc, window in zip(utcs, windows_channel) if len(window) > 0 and utc.year > 2020]

    # group time windows by hour
    dict_hour_window = defaultdict(list)
    for utc, windows in pairs_utc_window:
        dict_hour_window[(utc.date(), utc.hour)].append(windows)
    '''
    # print the result
    for (date, hour), windows in dict_hour_window.items():
        print(f'({date}, {hour}): {windows}')
    '''
    # create lists of dates with hours and number of transients/pulses
    date_hour = []
    num_pulses = []
    for (date, hour), windows in dict_hour_window.items():
        date_hour.append(datetime.combine(date, time(hour=hour)))
        windows = list(itertools.chain(*windows)) # merge all time windows in the same time unit
        num_pulses.append(len(windows))
    
    # plot the data
    fig, ax = plt.subplots()
    ax.plot(date_hour, num_pulses)

    # set the date format
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))

    # rotate date labels automatically
    fig.autofmt_xdate()

    plt.xlabel('Time (Date and Hour)')
    plt.ylabel('Number of Transients/Pulses')
    plt.title(f'Transient Rate Evolution of DU{du} in channel{channel}')
    plt.grid(True)
    plt.savefig(f'result/transient_rate_DU{du}_channel{channel}.png')

for du in du_list:
    # make dictionaries
    channels  = ['X', 'Y', 'Z']
    windows = {channel: {} for channel in channels}

    # read selected NPZ file
    npz_file = np.load(f'result/DU{du}_threshold5_separation100.npz', allow_pickle=True)
    gps_times = npz_file['gps_time']
    windows['X'] = npz_file['window_x']
    windows['Y'] = npz_file['window_y']
    windows['Z'] = npz_file['window_z']

    utcs = gps2utc(gps_times)

    for channel in channels:
        plot_transient_rate(utcs, windows, channel)

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