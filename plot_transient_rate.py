# import necessary packages
import argparse
from datetime import datetime, time
import glob
import grand.dataio.root_trees as rt # GRANDLIB
import numpy as np
np.set_printoptions(threshold=np.inf)
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
plt.rc('font', size=10)
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
import os

###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Search for transient in data traces and save information in NPZ files per DU.")

parser.add_argument('--path_to_data_dir',
                    dest='path_to_data_dir',
                    default='nancay_data/',
                    type=str,
                    help='Specifies the path of the directory containing the ROOT files to analyse.')

parser.add_argument("--file_specific",
                    dest="file_specific",
                    default='*.root',
                    type=str,
                    help="Specify which ROOT files to analyze.")

parse_args = parser.parse_args()

path_to_data_dir    = parse_args.path_to_data_dir
file_specific       = parse_args.file_specific

#####################################
# GET ROOT FILES AND DU INFORMATION #
#####################################

# check if the data directory or the file list is empty
if not os.path.exists(path_to_data_dir):
    raise FileNotFoundError(f'Data directory not found: {path_to_data_dir}')
root_files = sorted(glob.glob(os.path.join(path_to_data_dir, file_specific)))
if not root_files:
    raise ValueError("The provided file list is empty. Please ensure that you provide a non-empty list of ROOT files.")

# initiate TADC tree of the file list to get basic info
total_entries   = 0 # total number of entries across all files
max_dus         = 0 # maximum used DUs
id_max_dus_file = 0 # index of file that has maximum used DUs
for i, file in enumerate(root_files):
    tadc = rt.TADC(file)
    tadc.get_entry(0) # get the entry from current file
    total_entries += tadc.get_number_of_entries()
    cur_dus = len(tadc.get_list_of_all_used_dus()) # get used DUs from current file
    if cur_dus > max_dus:
        max_dus         = cur_dus
        id_max_dus_file = i

# get list of used DUs
tadc = rt.TADC(root_files[id_max_dus_file])
du_list = tadc.get_list_of_all_used_dus()
print('\nFiles contain data from following DUs: {}'.format(du_list))

###############################################
# COMPUTE AND PLOT TRANSIENT RATE FOR EACH DU #
###############################################

du = du_list[0]

# read the selected NPZ file
npz_file = np.load('result2/du_{}_threshold_5_separation_100.npz'.format(du), allow_pickle=True)
gps_time = npz_file['gps_time']
windows_x = npz_file['window_x']
windows_y = npz_file['window_y']
windows_z = npz_file['window_z']

# convert GPS time to UTC; GPS time = UTC + 18s at present
def gps2utc(gps_time):
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    utc_time = datetime.utcfromtimestamp(gps_time - leap_seconds)
    return utc_time
gps2utc_v = np.vectorize(gps2utc) # create a vectorized version of the converting function
utc_times = gps2utc_v(gps_time)

# compute and plot transient rate for 3 ADC channels of selected DU
def plot_transient_rate(utc_times, windows, channel):
    # use a list comprehension to filter out empty lists and pair time with non-empty time windows
    pairs_time_window = [(time, window) for time, window in zip(utc_times, windows) if len(window) > 0 and time.year > 2020]

    # group time windows by hour
    dict_hour_window = {}
    for utc_time, windows in pairs_time_window:
        if (utc_time.date(), utc_time.hour) not in dict_hour_window :
            dict_hour_window[(utc_time.date(), utc_time.hour)] = []
        dict_hour_window[(utc_time.date(), utc_time.hour)].append(windows)
    
    # print the result
    for (date, hour), windows in dict_hour_window.items():
        print(f'({date}, {hour}): {windows}')

    # create lists of dates with hours and number of transients/pulses
    date_hour = []
    num_pulses = []
    for (date, hour), windows in dict_hour_window.items():
        date_hour.append(datetime.combine(date, time(hour=hour)))
        windows = [window for entry in windows for window in entry] # merge time windows in the same time unit
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
    npz_file = 'du_{}_threshold_5_separation_100.npz'.format('{}')
    plt.title('Transient Rate Evolution of DU{} in channel{}'.format(du, channel))
    plt.grid(True)
    plt.savefig('result2/transient_rate_DU{}_channel{}.png'.format(du, channel))

plot_transient_rate(utc_times, windows_x, 'X')
plot_transient_rate(utc_times, windows_y, 'Y')
plot_transient_rate(utc_times, windows_z, 'Z')

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