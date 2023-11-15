#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#############################
# GET NPZ FILES AND DU LIST #
#############################

# get all NPZ files
file_list = sorted(glob.glob(npz_files))

# create an empty set to store unique DUs
du_set = set()

# show all NPZ files and extract DUs
print('\nData from following NPZ files:')
for file in file_list:
    basename = os.path.basename(file)
    print(basename)
    
    # split the filename on underscore
    du = basename.split('_')[0][2:]
    du_set.add(du)

# convert the set to a list in order
du_list = sorted(list(du_set))

# print the list of all used DUs
print(f'\nNPZ files contain data from following DUs: \n{du_list}\n')

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
    utc_hour_list = []
    rate_list     = []
    for (date, hour), windows in dict_hour_window.items():
        utc_hour_list.append(datetime.combine(date, time(hour=hour)))
        duration = len(windows) * num_samples * time_step # total time of traces
        windows  = list(itertools.chain(*windows)) # filter out empty time windows and reformat windows
        rate_list.append(len(windows) / duration * 1e6) # GHz --> kHz

    # convert UTCs to DunHuang local times
    DunHuang_hour_list = [utc_hour + timedelta(hours=8) for utc_hour in utc_hour_list]
    
    return DunHuang_hour_list, rate_list

# make dictionaries for easier indexing
hours_list = {channel: {} for channel in channels}
rates_list = {channel: {} for channel in channels}

# loop through all DUs to initiate dictionaries with empty lists
for du in du_list:
    for channel in channels:
        hours_list[channel][du] = []
        rates_list[channel][du] = []

# loop through all files to compute transient rates
for file in file_list:
    # get information from this file
    npz_file     = np.load(file, allow_pickle=True)
    du           = str(npz_file['du'])
    windows_list = {channel: npz_file[f'window{channel}'] for channel in channels}

    # convert GPS times to UTCs
    utcs = gps2utc(npz_file['gps_times'])

    # compute transient rates for 3 ADC channels
    for channel in channels:
        part_hours, part_rates = compute_rate(utcs, windows_list[channel])
        hours_list[channel][du].append(part_hours)
        rates_list[channel][du].append(part_rates)
    
# loop through all DUs to plot transient rates for 3 ADC channels
for du in du_list:
    # create a 3x1 layout for the subplots and adjust figure size
    fig, axes = plt.subplots(3, 1, figsize=(36, 12), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the data per channel
    for id, channel in enumerate(channels):
        # reformat DunHuang hours and transient rates to plot
        hours_list[channel][du] = list(itertools.chain(*hours_list[channel][du]))
        rates_list[channel][du] = 1+np.array(list(itertools.chain(*rates_list[channel][du])))

        # locate corresponding axis
        ax = axes[id]

        # scatter points and add vertical lines to make the plot more clearly
        ax.scatter(hours_list[channel][du], rates_list[channel][du], color='blue', s=9, alpha=0.6)
        ax.vlines(hours_list[channel][du], 0, rates_list[channel][du], color='blue', linestyle='solid', linewidth=4, alpha=0.75)

        # set axis scales
        ax.set_yscale('log')

        # set only the lower limit of the y-axis
        ax.set_ylim(bottom=1)

        # enable gird
        ax.grid(True)

        # enable ticks on both left and right sides of the plot
        ax.tick_params(axis='y', labelleft=True, labelright=True)

        # set title on the right-hand side
        ax.text(1.01, 0.5, f'Channel {channel}', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=20, rotation=-90)
        
        # set X-limits
        ax.set_xlim(datetime(2023, 10, 11, 6), datetime(2023, 10, 31, 18))

        # set the major tick locator to every 12 hours
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
        
        # add vertical lines at sunrise and sunset
        ylim = ax.get_ylim()
        ax.vlines([datetime(2023, 10, 11, 8, 00),
                   datetime(2023, 10, 12, 8, 00),
                   datetime(2023, 10, 13, 8, 00),
                   datetime(2023, 10, 14, 8, 00),
                   datetime(2023, 10, 15, 8, 00),
                   datetime(2023, 10, 16, 8, 00),
                   datetime(2023, 10, 17, 8, 00),
                   datetime(2023, 10, 18, 8, 00),
                   datetime(2023, 10, 19, 8, 00),
                   datetime(2023, 10, 20, 8, 00),
                   datetime(2023, 10, 27, 8, 00),
                   datetime(2023, 10, 28, 8, 00),
                   datetime(2023, 10, 29, 8, 00),
                   datetime(2023, 10, 30, 8, 00),
                   datetime(2023, 10, 31, 8, 00),], ymin=ylim[0], ymax=ylim[1], color='red', linestyles='dashed')
        ax.vlines([datetime(2023, 10, 11, 19, 00),
                   datetime(2023, 10, 12, 19, 00),
                   datetime(2023, 10, 13, 19, 00),
                   datetime(2023, 10, 14, 19, 00),
                   datetime(2023, 10, 15, 19, 00),
                   datetime(2023, 10, 16, 19, 00),
                   datetime(2023, 10, 17, 19, 00),
                   datetime(2023, 10, 18, 19, 00),
                   datetime(2023, 10, 19, 19, 00),
                   datetime(2023, 10, 20, 19, 00),
                   datetime(2023, 10, 26, 19, 00),
                   datetime(2023, 10, 27, 19, 00),
                   datetime(2023, 10, 28, 19, 00),
                   datetime(2023, 10, 29, 19, 00),
                   datetime(2023, 10, 30, 19, 00),
                   datetime(2023, 10, 31, 19, 00),], ymin=ylim[0], ymax=ylim[1], color='orange', linestyles='dashed')

        # only set X-tick labels for the last subplot
        if id < len(channels) - 1:
            ax.set_xticklabels([])
        else:
            # set the formatting of the ticks to include the hour
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))

            # rotate X labels
            for tick in ax.get_xticklabels():
                tick.set_rotation(15)

    # set X label and Y label for entire figure
    fig.text(0.5, 0.0, 'DunHuang Time in Date and Hour', ha='center', fontsize=20)
    fig.text(0.0, 0.5, 'Transient Rate / kHz', va='center', rotation='vertical', fontsize=20)
    
    # add the figure legend
    fig.legend(custom_lines, ['Sunrise', 'Sunset'], loc='upper right', fontsize=18, bbox_to_anchor=(0.98,1), bbox_transform=plt.gcf().transFigure)

    # add a main/super title for the entire figure
    plt.suptitle(f'Time Evolution of Transient Rates \n DU{du}, Threshold = {num_threshold}, Separation = {standard_separation}, Crossing = {num_crossings}, Cut-off frequency = {cutoff_frequency}, Noises = [{noises[0]}, {noises[1]}, {noises[2]}]', fontsize=24)

    # adjust layout
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 1.0])

    # save the figure as a PNG file
    rate_dir  = 'plot/rate/'
    rate_file = f'rate_DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_cutoff{cutoff_frequency}_noise{noises[0]}_{noises[1]}_{noises[2]}.png'
    plt.savefig(os.path.join(rate_dir, rate_file))
    print(f'Saved: {rate_dir}{rate_file}')

    # close the figure to free up memory
    plt.close(fig)

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")