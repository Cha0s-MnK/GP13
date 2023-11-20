#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

time_axis = np.arange(1024) * time_step # time axis of a trace

######################################
# PLOT TRACES TO CHECK SEARCH METHOD #
######################################

def plot_trace(file):
    # get information of this NPZ file
    npz_file = np.load(file, allow_pickle=True)
    traces   = {channel: npz_file[f'trace{channel}'] for channel in channels}
    utc      = npz_file['utc']
    utc      = datetime.strptime(utc.item(), '%Y%m%d%H%M%S')
    utc_flat = utc.strftime('%Y%m%d%H%M%S')

    # create a 6x1 layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(6, 1, figsize=(36, 24), dpi=200)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the original trace, the filtered trace and corresponding time windows per channel
    for id, channel in enumerate(channels):
        # original traces will be on the first three subplots: 0, 1, 2
        plot_channel(axes[id], traces, channel, 'off')

        # filtered traces will be on the last three subplots: 3, 4, 5
        plot_channel(axes[id+3], traces, channel, 'on')

    # set common labels and title
    fig.text(0.5, 0.0, 'Time / ns', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'ADC Counts', va='center', rotation='vertical', fontsize=18)
    plt.suptitle(f'Original Traces and Filtered Traces with Corresponding Searched Windows\n{utc}', fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])

    # save the figure as a PNG file
    check_plot_name = f'check_DU{good_du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{utc_flat}.png'
    check_plot_file = os.path.join(check_plot_dir, check_plot_name)
    plt.savefig(check_plot_file)
    print(f'Saved: {check_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

def plot_channel(ax, traces, channel, filter_status):
    # search for time windows
    threshold   = num_threshold * noises[channel][good_du]
    window_list = search_windows(trace=traces[channel], 
                                 threshold=threshold, 
                                 filter=filter_status)

    # filter trace based on filter status
    if filter_status == 'on':
        trace = high_pass_filter(traces[channel], sample_frequency, cutoff_frequency)
    else:
        trace = traces[channel]

    # plot the trace
    ax.plot(time_axis, trace, color='blue', label=f'trace')

    # add horizontal dashed lines at +/- threshold of the trace
    ax.axhline(y=threshold, color='red', linestyle='--', label=f'+/- threshold={threshold}')
    ax.axhline(y=-threshold, color='red', linestyle='--')

    # add horizontal dashed lines at +/- STD of the trace
    std = np.std(trace)
    ax.axhline(y=std, color='orange', linestyle='--', label=f'+/- STD={std:.2f}')
    ax.axhline(y=-std, color='orange', linestyle='--')


    # highlight the time windows
    for id, window in enumerate(window_list):
        start_id, stop_id = window
        start_time        = start_id * time_step
        stop_time         = stop_id * time_step
        if id == 0:
            ax.axvspan(start_time, stop_time, color='green', alpha=0.3, label='Time Window(s)')
        else:
            ax.axvspan(start_time, stop_time, color='green', alpha=0.3)

        # annotate time windows
        centre_time = (start_time + stop_time) / 2
        ylim = ax.get_ylim()
        ax.text(centre_time, 1.02*ylim[1], f'[{start_time}, {stop_time}]', color='black', fontsize=10, ha='center', va='bottom')

    # set title on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel}, filter {filter_status}', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    # set X-limits
    #ax.set_xlim(-4, 4100)

    # set X-tick every 250 units
    ax.set_xticks(np.arange(min(time_axis), max(time_axis)+1, 250))

    # add legends
    ax.legend(loc='upper right')

    # enable grid
    ax.grid(True)

    # enable Y-ticks on both sides
    ax.tick_params(axis='y', labelleft=True, labelright=True)

#################
# MAIN FUNCTION #
#################

def main():
    # get NPZ files
    file_list = sorted(glob(check_plot_files))

    # check 1 ROOT file to test the search method; only consider 1 good DU
    for file in file_list:
        plot_trace(file)
    pass

if __name__ == "__main__":
    with record_run_time():
        main()