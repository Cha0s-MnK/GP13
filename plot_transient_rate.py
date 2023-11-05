# import necessary packages
import numpy as np
np.set_printoptions(threshold=np.inf)
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

###############################################
# COMPUTE AND PLOT TRANSIENT RATE FOR EACH DU #
###############################################

# read the selected NPZ file
npz_file = np.load('result/du_1013_threshold_5_separation_100.npz', allow_pickle=True)
gps_time = npz_file['gps_time']
windows_x = npz_file['window_x']
windows_y = npz_file['window_y']
windows_z = npz_file['window_z']

# convert GPS time to UTC; GPS time = UTC + 18s at present
from datetime import datetime, time
def gps2utc(gps_time):
    leap_seconds = 18 # number of leap seconds since Jan 6th 1980
    utc_time = datetime.utcfromtimestamp(gps_time - leap_seconds)
    return utc_time
gps2utc_v = np.vectorize(gps2utc) # create a vectorized version of the converting function
utc_times = gps2utc_v(gps_time)

# compute and plot transient rate for 3 ADC channels of selected DU
def plot_transient_rate(utc_times, windows):
    # use a list comprehension to filter out empty lists and pair time with non-empty time windows
    pairs_time_window = [(time, window) for time, window in zip(utc_times, windows) if len(window) > 0 and time.year > 2020]

    # group time windows by hour
    dict_hour_window = {}
    for time, window in pairs_time_window:
        if (time.date(), time.hour) not in dict_hour_window :
            dict_hour_window[(time.date(), time.hour)] = []
        dict_hour_window[(time.date(), time.hour)].append(window)
    '''
    # merge windows in the same hour
    for (date, hour), windows in dict_hour_window.items():
        windows = [window for sublist in windows for window in sublist] # flatten time windows
        dict_hour_window[(date, hour)] = windows
    '''
    # print the result
    for (date, hour), windows in dict_hour_window.items():
        print(f'({date}, {hour}): {windows}')

plot_transient_rate(utc_times, windows_x)
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