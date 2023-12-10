#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

# create the save directory if it does not exist
name = 'rate'
save_dir = os.path.join(plot_dir, name)
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

def plot1du_rate(du):
    file_list = get_npz_files(name=name, date_list=False, du_list=[du])

    # make dictionaries for easier indexing and initiate them
    hours_list        = {channel: [] for channel in channels}
    rates_list        = {channel: [] for channel in channels}
    filter_hours_list = {channel: [] for channel in channels}
    filter_rates_list = {channel: [] for channel in channels}

    for file in file_list:
        # get info from NPZ files
        npz_file            = np.load(file, allow_pickle=True)
        windows_list        = {channel: npz_file[f'windows_list{channel}'] for channel in channels}
        filter_windows_list = {channel: npz_file[f'filter_windows_list{channel}'] for channel in channels}
        csts_list           = {channel: npz_file[f'csts_list{channel}'] for channel in channels}

        # compute transient/pulse rates
        for channel in channels:
            hour_list, rate_list = get1channel_rate(csts_list[channel], windows_list[channel])
            hours_list[channel].append(hour_list)
            rates_list[channel].append(rate_list)
            hour_list, rate_list = get1channel_rate(csts_list[channel], filter_windows_list[channel])
            filter_hours_list[channel].append(hour_list)
            filter_rates_list[channel].append(rate_list)

    # reformat DunHuang hours and transient rates
    for channel in channels:
        hours_list[channel]        = list(itertools.chain(*hours_list[channel]))
        rates_list[channel]        = list(itertools.chain(*rates_list[channel]))
        rates_list[channel]        = [rate + 0.1 for rate in rates_list[channel]] # add 0.1Hz to distinguish it from no data taken

        filter_hours_list[channel] = list(itertools.chain(*filter_hours_list[channel]))
        filter_rates_list[channel] = list(itertools.chain(*filter_rates_list[channel]))
        filter_rates_list[channel] = [rate + 0.1 for rate in filter_rates_list[channel]]

    fig, axes = plt.subplots(3, 1, figsize=(24, 12), dpi=165)

    # plot the transient/pulse rates of the filtered trace per channel
    for id, channel in enumerate(channels):
        plot1channel_Dec(ax=axes[id], channel=channel, hour_list=hours_list[channel], rate_list=rates_list[channel])

    # set common labels and title
    fig.text(0.5, 0.01, 'CST', ha='center', fontsize=18)
    fig.text(0.01, 0.5, 'Transient/Pulse Rate / Hz', va='center', rotation='vertical', fontsize=18)
    #plt.suptitle(f'Time Evolution of Transient Rates for DU{du}', fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.02, 0.02, 0.95, 1.0])

    # save the figure as a PNG file
    save_file = os.path.join(save_dir, f'rate_DU{du}_RUN{nums_run}_T1st{num_threshold1}_T2nd{num_threshold2}_MIN{num_crossing_min}_MAX{num_crossing_max}_interval{num_interval}cutoff{cutoff_frequency}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')
    plt.close(fig)

# compute DunHuang hours and transient/pulse rates
def get1channel_rate(cst_list, window_list):
    # pair UTCs with transient/pulse windows
    pairs_cst_window = [(cst, window) for cst, window in zip(cst_list, window_list)]

    # make a dictionary to pair merge time windows with the same hour together
    dict_hour_windows = defaultdict(list)
    for cst, window in pairs_cst_window:
        dict_hour_windows[(cst.date(), cst.hour)].append(window)

    # build lists of UTC hours and transient rates
    hour_list = []
    rate_list = []
    for (date, hour), windows in dict_hour_windows.items():
        hour_list.append(datetime.combine(date, time(hour=hour)))
        duration = len(windows) * num_samples * time_step # total time of traces
        windows  = list(itertools.chain(*windows)) # filter out empty time windows and reformat windows
        rate_list.append(len(windows) / duration * 1e9) # [GHz] --> [Hz]

    return hour_list, rate_list

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

def plot1channel_Dec(ax, channel, hour_list, rate_list):
    print(np.mean(rate_list))
    # scatter points and add vertical lines to make the plot more clear
    ax.scatter(hour_list, rate_list, color='blue', s=9, alpha=0.6)
    ax.vlines(hour_list, 0, rate_list, color='blue', linestyle='solid', linewidth=4, alpha=0.75)

    # set Y-scale
    ax.set_yscale('log')
    #ax.set_ylim([5e-2, 1e6])

    # set the title on the right-hand side
    ax.text(1.01, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)
    
    # set X-ticks
    ax.set_xlim(datetime(2023, 12, 5, 12), datetime(2023, 12, 8, 12))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    for tick in ax.get_xticklabels():
        tick.set_rotation(30)
    
    # add horizontal line at 678Hz (1 transient/pulse in 1 hour)
    xlim = ax.get_xlim()
    ax.hlines(678, xmin=xlim[0], xmax=xlim[1], color='green', linestyle=':', label='678Hz')

    # add vertical lines at sunrise and sunset
    ylim = ax.get_ylim()
    ax.vlines([datetime(2023, 12, 6, 8, 30),
               datetime(2023, 12, 7, 8, 30),
               datetime(2023, 12, 8, 8, 30)], 
               ymin=ylim[0], 
               ymax=ylim[1],
               color='red', 
               linestyle='dashed', 
               label='sunrise')
    ax.vlines([datetime(2023, 12, 5, 18, 30),
               datetime(2023, 12, 6, 18, 30),
               datetime(2023, 12, 7, 18, 30),
               datetime(2023, 12, 8, 18, 30)], 
               ymin=ylim[0], 
               ymax=ylim[1], 
               color='orange', 
               linestyles='dashed', 
               label='sunset')

    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', fontsize=16)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

#################
# MAIN FUNCTION #
#################

def main():
    du_list = get_npz_dus(file_list=get_npz_files(name=name, du_list=False))

    for du in du_list:
        plot1du_rate(du=du)

if __name__ == "__main__":
    with record_run_time():
        main()