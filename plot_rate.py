#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

################################################
# COMPUTE AND PLOT TRANSIENT RATES FOR EACH DU #
################################################

def plot_du_rate():
    # make dictionaries for easier indexing and initiate them with empty lists
    hours_list          = {channel: [] for channel in channels}
    rates_list          = {channel: [] for channel in channels}
    filtered_hours_list = {channel: [] for channel in channels}
    filtered_rates_list = {channel: [] for channel in channels}

    # loop through files to get DunHuang hours and transient rates
    for file in file_list:
        # check DU
        basename = os.path.basename(file)
        if int(basename.split('_')[1][2:]) != du:
            continue

        # get information from this file
        npz_file              = np.load(file, allow_pickle=True)
        windows_list          = {channel: npz_file[f'windows_list{channel}'] for channel in channels}
        filtered_windows_list = {channel: npz_file[f'filtered_windows_list{channel}'] for channel in channels}
        gps_time_list         = npz_file['gps_time_list']

        # convert GPS times to UTCs
        utc_list = gps2utc(gps_time_list)

        # compute transient/pulse rates for 3 ADC channels
        for channel in channels:
            hour_list, rate_list = compute_rate(utc_list, windows_list[channel])
            hours_list[channel].append(hour_list)
            rates_list[channel].append(rate_list)
            hour_list, rate_list = compute_rate(utc_list, filtered_windows_list[channel])
            filtered_hours_list[channel].append(hour_list)
            filtered_rates_list[channel].append(rate_list)

        # reformat DunHuang hours and transient rates
        for channel in channels:
            hours_list[channel]          = list(itertools.chain(*hours_list[channel]))
            rates_list[channel]          = list(itertools.chain(*rates_list[channel]))
            # add 0.1 to distinguish it from no data taken
            rates_list[channel]          = [rate + 0.1 for rate in rates_list[channel]]
            filtered_hours_list[channel] = list(itertools.chain(*filtered_hours_list[channel]))
            filtered_rates_list[channel] = list(itertools.chain(*filtered_rates_list[channel]))
            filtered_rates_list[channel] = [rate + 0.1 for rate in filtered_rates_list[channel]]

        # plot transient/pulse rates
        plot_rate(du, hours_list, rates_list, filtered_hours_list, filtered_rates_list)

def get_du_rate_Nov(du):
    # update NPZ files according to current DU
    file_list = sorted(glob(os.path.join(save_dir, f'*/*DU{du}_RUN{num_run}*.npz')))

    # make dictionaries for easier indexing and initiate them with empty lists
    hours_list = {channel: [] for channel in channels}
    rates_list = {channel: [] for channel in channels}

    # loop through NPZ files to get DunHuang hours and transient rates
    for file in file_list:
        # get information from this file
        npz_file      = np.load(file, allow_pickle=True)
        windows_list  = {channel: npz_file[f'filtered_windows_list{channel}'] for channel in channels}
        gps_time_list = npz_file['gps_time_list']

        # convert GPS times to UTCs
        utc_list = gps2utc(gps_time_list)

        # compute transient/pulse rates for 3 ADC channels
        for channel in channels:
            hour_list, rate_list = compute_rate(utc_list, windows_list[channel])
            hours_list[channel].append(hour_list)
            rates_list[channel].append(rate_list)

    # reformat DunHuang hours and transient rates
    for channel in channels:
        hours_list[channel] = list(itertools.chain(*hours_list[channel]))
        rates_list[channel] = list(itertools.chain(*rates_list[channel]))

        # add 0.1 to distinguish it from no data taken
        rates_list[channel] = [rate + 0.1 for rate in rates_list[channel]]

    # plot transient/pulse rates
    plot_du_rate_Nov(du, hours_list, rates_list)
    
    pass

# compute DunHuang hours and transient/pulse rates
def compute_rate(utc_list, windows_channel):
    # pair UTCs with time windows
    pairs_utc_window = [(utc, window) for utc, window in zip(utc_list, windows_channel)]

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
        # [GHz] --> [Hz]
        rate_list.append(len(windows) / duration * 1e9)

    # convert UTCs to DunHuang local times
    DunHuang_hour_list = [utc_hour + timedelta(hours=8) for utc_hour in utc_hour_list]

    return DunHuang_hour_list, rate_list

def plot_rate(du, hours_list, rates_list, filtered_hours_list, filtered_rates_list):
    # create a 6x1 layout for the subplots and adjust the figure size
    fig = plt.figure(figsize=(36, 24), dpi=256)
    gs = gridspec(6, 1, figure=fig)

    # plot the transient/pulse ratesof the original trace and the filtered trace per channel
    for id, channel in enumerate(channels):
        # original traces will be on the first three subplots: 0, 1, 2
        plot_channel(gs, id, hours_list, rates_list, channel, 'off')

        # filtered traces will be on the last three subplots: 3, 4, 5
        plot_channel(gs, id+3, filtered_hours_list, filtered_rates_list, channel, 'on')

    # set common labels and title
    fig.text(0.5, 0.0, 'DunHuang Time in Date and Hour', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'Transient or Pulse Rate / kHz', va='center', rotation='vertical', fontsize=18)
    plt.suptitle(f'Time Evolution of Transient Rates for DU{du}, \nThreshold = {num_threshold}, Max Separation = {standard_separation}, Min Crossing = {num_crossings}, Max Samples = {max_samples}, Max Fluctuation = {fluctuation}, Sample Frequency = {sample_frequency}, Cutoff frequency = {cutoff_frequency}', fontsize=20)

    # add legend
    fig.legend(custom_sun, ['Sunrise', 'Sunset'], loc='upper right', fontsize=18, bbox_to_anchor=(0.98,1), bbox_transform=plt.gcf().transFigure)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])

    # save the figure as a PNG file
    rate_plot_name = f'rate_DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_fluctuation{fluctuation}_frequency{sample_frequency}_cutoff{cutoff_frequency}.png'
    rate_plot_file = os.path.join(rate_plot_dir, rate_plot_name)
    plt.savefig(rate_plot_file)
    print(f'Saved: {rate_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

def plot_du_rate_Nov(du, hours_list, rates_list):
    # create a 3x1 layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(3, 1, figsize=(36, 12), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the transient/pulse rates of the filtered trace per channel
    for id, channel in enumerate(channels):
        plot_channel_Nov(axes[id], hours_list[channel], rates_list[channel], channel)

    # set common labels and title
    fig.text(0.5, 0.0, 'DunHuang Time in Date and Hour', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'Transient or Pulse Rate / Hz', va='center', rotation='vertical', fontsize=18)
    plt.suptitle(f'Time Evolution of Transient Rates for DU{du}, \nThreshold = {num_threshold}, Max Separation = {standard_separation}, Min Crossing = {num_crossings}, Max Samples = {max_samples}, Max Fluctuation = {std_fluctuation}, Sample Frequency = {sample_frequency}, Cutoff frequency = {cutoff_frequency}', fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])

    # save the figure as a PNG file
    rate_plot_file = os.path.join(rate_plot_dir, f'rate_DU{du}_RUN{num_run}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_fluctuation{std_fluctuation}_frequency{sample_frequency}_cutoff{cutoff_frequency}.png')
    plt.savefig(rate_plot_file)
    print(f'Saved: {rate_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

def plot_channel(ax, 
                 hours_list, 
                 rates_list, 
                 channel, 
                 filter_status):
    # scatter points and add vertical lines to make the plot more clearly
    bax.scatter(hours_list[channel], rates_list[channel], color='blue', s=9, alpha=0.6)
    bax.vlines(hours_list[channel], 0, rates_list[channel], color='blue', linestyle='solid', linewidth=4, alpha=0.75)

    # enable gird
    bax.grid(True)

    # set Y-scale
    bax.set_yscale('log')

    # enable Y-ticks on both sides
    bax.tick_params(axis='y', labelleft=True, labelright=True)

    # set the title on the right-hand side
    bax.text(1.01, 0.5, f'channel {channel}, filter {filter_status}', verticalalignment='center', horizontalalignment='left', transform=bax.transAxes, fontsize=16, rotation=-90)
    
    # set X-ticks
    bax.set_xlim(datetime(2023, 11, 19, 0), datetime(2023, 11, 22, 0))
    bax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    bax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    for tick in bax.get_xticklabels():
        tick.set_rotation(15)
        
    # add vertical lines at sunrise and sunset
    ylim = bax.get_ylim()
    bax.vlines([datetime(2023, 10, 11, 8, 00),
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
    bax.vlines([datetime(2023, 10, 11, 19, 00),
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

def plot_channel_Nov(ax, 
                     hour_list, 
                     rate_list, 
                     channel):
    # scatter points and add vertical lines to make the plot more clear
    ax.scatter(hour_list, rate_list, color='blue', s=9, alpha=0.6)
    ax.vlines(hour_list, 0, rate_list, color='blue', linestyle='solid', linewidth=4, alpha=0.75)

    # set Y-scale
    ax.set_yscale('log')
    ax.set_ylim([5e-2, 1e6])

    # set the title on the right-hand side
    ax.text(1.01, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)
    
    # set X-ticks
    ax.set_xlim(datetime(2023, 11, 19, 0), datetime(2023, 11, 25, 0))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    for tick in ax.get_xticklabels():
        tick.set_rotation(30)
    
    # add horizontal line at 678Hz (1 transient/pulse in 1 hour)
    xlim = ax.get_xlim()
    ax.hlines(678, xmin=xlim[0], xmax=xlim[1], color='green', linestyle=':', label='1 in 1 hour')

    # add vertical lines at sunrise and sunset
    ylim = ax.get_ylim()
    ax.vlines([datetime(2023, 11, 16, 8, 30),
               datetime(2023, 11, 17, 8, 30),
               datetime(2023, 11, 18, 8, 30),
               datetime(2023, 11, 19, 8, 30),
               datetime(2023, 11, 20, 8, 30),
               datetime(2023, 11, 21, 8, 30),
               datetime(2023, 11, 22, 8, 30),
               datetime(2023, 11, 23, 8, 30),
               datetime(2023, 11, 24, 8, 30)], 
               ymin=ylim[0], 
               ymax=ylim[1], 
               color='red', 
               linestyle='dashed', 
               label='sunrise')
    ax.vlines([datetime(2023, 11, 15, 18, 30),
               datetime(2023, 11, 16, 18, 30),
               datetime(2023, 11, 17, 18, 30),
               datetime(2023, 11, 18, 18, 30),
               datetime(2023, 11, 19, 18, 30),
               datetime(2023, 11, 20, 18, 30),
               datetime(2023, 11, 21, 18, 30),
               datetime(2023, 11, 22, 18, 30),
               datetime(2023, 11, 23, 18, 30),
               datetime(2023, 11, 24, 18, 30)], 
               ymin=ylim[0], 
               ymax=ylim[1], 
               color='orange', 
               linestyles='dashed', 
               label='sunset')

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, loc='upper fight', fontsize=16)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

    pass


#################
# MAIN FUNCTION #
#################

def main():
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(save_dir, f'*/*RUN{num_run}*.npz')))

    # get used DUs
    du_list = get_npz_dus(file_list)

    # loop through used DUs to compute and plot transient/pulse rates
    for du in du_list:
        get_du_rate_Nov(du)

if __name__ == "__main__":
    with record_run_time():
        main()