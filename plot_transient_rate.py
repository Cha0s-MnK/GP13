# import necessary packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
plt.rc('font', size=10)
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# system
import os
import glob
import argparse
import datetime

# scipy
import numpy as np

# grandlib
import grand.dataio.root_trees as rt

################################################
# SEARCHING FOR TRANSIENTS IN GRAND DATA FILES #
################################################

# GRANDLIB BRANCH: 'dev_io_root'
# This script reads out a GRAND data file in ROOT format. It searches for transient pulses in the time traces of the data 
# file, where each DU is treated separately. Here, a transient is defined as a pulse that exceeds 5 times the STD of the 
# trace, called 5 sigma from this point onwards. For each trace, the number of 5 sigma transients is computed, as well as 
# the number of times this 5 sigma threshold is exceeded during the transient. The maximum peak position of the transient is
# also stored. All information is saved in an NPZ file for each DU.

###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Search for transient in data traces and save information in NPZ files per DU.")

parser.add_argument('--path_to_data_dir',
                    dest='path_to_data_dir',
                    default='gp13_data/',
                    type=str,
                    help='Specifies the path of the directory containing the ROOT files to analyse.')

parser.add_argument("--file_specific",
                    dest="file_specific",
                    default='*.root',
                    type=str,
                    help="Specify which ROOT files to analyze.")

parser.add_argument('--empty_channel',
                    dest='empty_channel',
                    default=0,
                    type=int,
                    help='Specifies which of the 4 ADC channels is not branched.')

parser.add_argument('--num_threshold',
                    dest='num_threshold',
                    default=5,
                    type=int,
                    help='Specifies the trigger threshold in units of sigma/STDs.')

parser.add_argument('--standard_separation',
                    dest='standard_separation',
                    default=100,
                    type=int,
                    help='Specifies the required separation between two pulses in a trace in units of sample numbers.')                    

parse_args = parser.parse_args()

path_to_data_dir    = parse_args.path_to_data_dir
file_specific       = parse_args.file_specific
empty_channel       = parse_args.empty_channel
num_threshold       = parse_args.num_threshold
standard_separation = parse_args.standard_separation

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
'''
# Sometimes all four channels are enabled, but one is not branched. For GP13 it is always channel 0
mask_channel = np.array(tadc.adc_enabled_channels_ch[0], dtype=bool)
if np.all(mask_channel):
    mask_channel[empty_channel] = False

# Sometimes there seems to be a bug where there are less than 3 enabled ADC channels.
if len(mask_channel[mask_channel]) < 3:
    mask_channel = np.ones(4, dtype=bool)
    mask_channel[empty_channel] = False
'''
# create a channel mask to enable 3 ADC channels and disable channel 0 for GP13
mask_channel = np.ones(4, dtype=bool)
mask_channel[empty_channel] = False
'''
# get some useful parameters from the data directory
path_to_data_dir_split = path_to_data_dir.split('/')
array = path_to_data_dir_split[4]
month = path_to_data_dir_split[5]
mode  = path_to_data_dir_split[6]

# define the NPZ directory and file name template to store the computed PSDs
npz_dir           = 'data/{}/{}/{}/transient_search/'.format(array, month, mode)
npz_file_template = 'transient_search_{}_{}_{}_du_{}_thresh_{}_separ_{}.npz'.format(array, month, mode, '{}', thresh, standard_separation)
'''
# define the NPZ directory and file name template to store the computed PSDs
npz_dir           = 'result/'
npz_file_template = 'du_{}_threshold_{}_separation_{}.npz'.format('{}', num_threshold, standard_separation)

def search_transients(trace,
                      num_threshold=5,
                      standard_separation=100): 
    '''
    Identify transients of a signal where the signal exceeds a certain threshold.

    Parameters:
    trace (array-like): the signal to be analyzed
    num_threshold (float): a multiplier for the STD to determine the threshold for identifying transients.
    standard_separation (int): a minimum separation between transients in units of samples.

    Returns:
    windows: a list of windows, each containing the start and end points of a transient.
    crossings: a list of crossings, each containing the values of the signal during a transient.
    '''
    # create the return lists
    windows       = [] # time windows
    crossings     = [] # threshold crossings
    traces_window = [] # traces in time windows

    # find trace positions where N sigma/STD is exceeded
    threshold         = num_threshold * np.std(trace)
    exceed_threshold  = np.abs(trace) > threshold
    indices_crossings = np.nonzero(exceed_threshold)[0]

    # stop if there are no transients
    if not np.any(exceed_threshold):
        return windows, crossings, traces_window

    # find the separations between consecutive threshold crossings
    separations_crossings = np.diff(indices_crossings)

    # locate where two transients/pulses are separated
    indices_pulses = np.concatenate((np.nonzero(separations_crossings > standard_separation)[0], [len(indices_crossings)-1]))
    
    # search all transients/pulses
    half_separation = standard_separation / 2
    for index in indices_pulses:
        # get the beginning of current pulse
        if index == indices_pulses[0]:
            begin_pulse = np.max([0, indices_crossings[0] - half_separation])
        else:
            begin_pulse = end - standard_separation + separations_crossings[index-1]

        # get the end of current pulse
        end = np.min([len(trace)-1, indices_crossings[index] + half_separation] )

        # count how many crossings (in units of sigma/STD) are in the time window
        abs_trace_window = np.abs(trace[int(begin_pulse) : int(end)+1])
        crossing         = [signal * num_threshold / threshold for signal in abs_trace_window if signal > threshold]

        # add info of current pulse to the return lists
        windows.append([begin_pulse, end])
        crossings.append(list(crossing))
        traces_window.append(trace[int(begin_pulse) : int(end)+1])

    return windows, crossings, traces_window

#########################################
# COMPUTE AND PLOT MEAN PSD FOR EACH DU #
#########################################

# make dictionaries so that you don't have to loop over all files each time you want to analyze a different DU
channels  = ['x', 'y', 'z']
traces = {channel: {} for channel in channels}
windows = {channel: {} for channel in channels}
crossings = {channel: {} for channel in channels}
traces_windows = {channel: {} for channel in channels}
gps_time = {}
for du in du_list:
    for channel in channels:
        windows[channel][du]   = np.zeros((total_entries), dtype=object) # save time windows of traces (unit = sample number)
        crossings[channel][du] = np.zeros((total_entries), dtype=object) # save threshold crossings of transients (unit = sigma/STD)
        traces_windows[channel][du]   = np.zeros((total_entries), dtype=object)
    gps_time[du] = np.zeros((total_entries), dtype=int) # save GPS times

# loop over all files in run
num_files = len(root_files)
id_entry  = 0
for i, root_file in enumerate(root_files):
    tadc  = rt.TADC(root_file)
    trawv = rt.TRawVoltage(root_file)
    num_entries = tadc.get_number_of_entries()

    # loop over all events in file
    print('\n{}/{}: Looping over {} events in {}'.format(i+1, num_files, num_entries, root_file))
    for entry in range(num_entries):
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # use the DU as key for the dictionaries
        # Assume that there are no coincident events, i.e. only one DU (and trace etc) per event.
        du = tadc.du_id[0]
    
        # get traces
        new_traces = np.array(tadc.trace_ch[0], dtype=object)[mask_channel]
        for channel, new_trace in zip(channels, new_traces):
            traces[channel] = new_trace
        '''
        # check if sample numbers are enough
        if any(len(channel) != num_samples for channel in new_traces[:3]):
            print('\nWARNING for entry {}: {} samples in trace != (3,{})'.format(entry, new_traces.shape, num_samples) + '\nSkipping...')
            id_entry += 1
            continue
        
        # check if data is between 10-55ºC
        if trawv.gps_temp[0] > 55 or trawv.gps_temp[0] < 10:
            print('\nWARNING for entry {}: gps_temp = {} not between 10 and 55 ºC'.format(entry,trawv.gps_temp[0]) + '\nSkipping...')
            id_entry += 1
            continue
        '''
        for channel in channels:
            new_windows, new_crossings, new_traces_windows = search_transients(traces[channel],
                                                                               num_threshold=num_threshold,
                                                                               standard_separation=standard_separation)
            windows[channel][du][id_entry] = new_windows
            crossings[channel][du][id_entry] = new_crossings
            traces_windows[channel][du][id_entry] = new_traces_windows
        gps_time[du][id_entry] = trawv.gps_time[0]

        id_entry += 1

for du in du_list:
    # get rid of events that are not filled because of bad quality data
    mask_filled = np.where(windows['x'][du] != 0)
    for channel in channels:
        windows[channel][du] = windows[channel][du][mask_filled]
        crossings[channel][du] = crossings[channel][du][mask_filled]
        traces_windows[channel][du] = traces_windows[channel][du][mask_filled]
    gps_time[du] = gps_time[du][mask_filled]

    # save total power in NPZ file
    npz_file = npz_file_template.format(du)
    np.savez(f'{npz_dir}{npz_file}',
             du = du,
             gps_time = gps_time[du],
             window_x = windows['x'][du],
             window_y = windows['y'][du],
             window_z = windows['z'][du],
             crossing_x = crossings['x'][du],
             crossing_y = crossings['y'][du],
             crossing_z = crossings['z'][du],
             trace_window_x = traces_windows['x'][du],
             trace_window_y = traces_windows['y'][du],
             trace_window_z = traces_windows['z'][du])
    print(f'\nSaved NPZ file: {npz_dir}{npz_file}')

# read the selected NPZ file
npz_file = np.load('result/du_1013_threshold_5_separation_100.npz', allow_pickle=True)
gps_time = npz_file['gps_time']
windows_x = npz_file['window_x']
'''
# use a list comprehension to filter out empty lists and flatten it
window_x = np.array([a_list for a_list in npz_file['window_x'] if len(a_list) > 0])
window_x_flat = [pair for sublist in window_x for pair in sublist]

# sort the array by the 1st element of each time window
window_x_sort = sorted(window_x_flat, key=lambda x: x[0])

# Print the sorted array
print(window_x_sort)
for pair in window_x_sort:
    print(pair)
'''
# convert GPS time to UTC; GPS time = UTC + 18s at present
from datetime import datetime, time

def gps2utc(gps_time):
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    utc_time = datetime.utcfromtimestamp(gps_time - leap_seconds)
    return utc_time

# create a vectorized version of the converting function
gps2utc_v = np.vectorize(gps2utc)

# convert all GPS times to UTC
utc_times = gps2utc_v(gps_time)

# pair UTC hours with time window and sort them
pairs_hour_window = [(time.date(), time.hour, window) for time, window in zip(utc_times, windows_x)]
pairs_hour_window.sort()

# group time windows by hour
dict_hour_window = {}
for date, hour, window in pairs_hour_window:
    if (date, hour) not in dict_hour_window :
        dict_hour_window[(date, hour)] = []
    dict_hour_window[(date, hour)].append(window)

# merge windows in the same hour
for (date, hour), windows in dict_hour_window.items():
    windows = np.array([window for window in windows if len(window) > 0]) # filter out empty time windows
    windows = [window for sublist in windows for window in sublist] # flatten time windows
    #windows = sorted(windows, key=lambda x: x[0]) # sort time windows by the 1st element
    dict_hour_window[(date, hour)] = windows

# remove entries with empty time windows
dict_hour_window = {key: value for key, value in dict_hour_window.items() if len(value) > 0}
'''
# print the result
for (date, hour), windows in dict_hour_window.items():
    print(f'({date}, {hour}): {windows}')
'''
# create lists of dates with hours and number of windows
date_hour = []
num_pulses = []
for (date, hour), windows in dict_hour_window.items():
    if date.year > 2020: # fix data bugs
        date_hour.append(datetime.combine(date, time(hour=hour)))
        num_pulses.append(len(windows))

print(date_hour)
print(num_pulses)
'''
# plot the data
fig, ax = plt.subplots()
ax.plot(date_hour, num_pulses)

# set the date format
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))

# rotate date labels automatically
fig.autofmt_xdate()

plt.xlabel('Time (Date and Hour)')
plt.ylabel('Number of Transients/Pulses')
plt.title('Evolution of Transient Rate')
plt.grid(True)
plt.savefig('result/test2.png')
'''
'''
# convert UTC to UTC hours and pair them with time windows
hour_window_pairs = [(datetime.utcfromtimestamp(time).hour, window) for time, window in zip(utc_times, window_x)]
print(hour_window_pairs)

# count the occurrences of each pair
counter = Counter(hour_window_pairs)

# Print the counter
for pair, count in counter.items():
    print(f"{pair}: {count}")
'''
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