#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

######################################
# PLOT TRACES TO CHECK SEARCH METHOD #
######################################

def check(file, good_du):
    # get info of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'Loop through the file with {num_entries} entries: \n{file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # 1 entry corresponds to only 1 DU; only consider the good DU
        if tadc.du_id[0] != good_du:
            continue
        
        # get the traces and plot
        traces = np.array(tadc.trace_ch[0])[channel_mask]
        gps_time = trawv.gps_time[0]
        plot_traces(traces, gps_time, good_du)
    pass

def plot_traces(traces, gps_time, good_du):
    # convert GPS time to UTC
    utc_time = gps2utc1(gps_time)

    # create a 3x1 layout for the subplots and adjust figure size
    fig, axes = plt.subplots(6, 1, figsize=(36, 24), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the original trace, the filtered trace and corresponding time windows per channel
    for id, channel in enumerate(channels):
        # original traces will be on the first three subplots: 0, 1, 2
        plot_channel(axes[id], traces[id], channel, 'off', good_du)

        # filtered traces will be on the last three subplots: 3, 4, 5
        plot_channel(axes[id+3], traces[id], channel, 'on', good_du)

    # set common labels and title
    fig.text(0.5, 0.0, 'Time / ns', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'ADC Counts', va='center', rotation='vertical', fontsize=18)
    plt.suptitle(f'Original Traces and Filtered Traces with Corresponding Searched Windows\n{utc_time}', fontsize=20)

    # adjust layout
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 1.0])

    # save the figure as a PNG file
    check_plot_name = f'DU{good_du}_{utc_time}.png'
    check_plot_file = os.path.join(check_plot_dir, check_plot_name)
    plt.savefig(check_plot_file)
    print(f'Saved: {check_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

def plot_channel(ax, trace, channel, filter_status, good_du):
    # filter trace based on filter status
    if filter_status == 'on':
        trace = high_pass_filter(trace, sample_frequency, cutoff_frequency)

    # search for time windows
    threshold = num_threshold*noises[channel][good_du]
    window_list = search_windows(trace=trace, 
                                 threshold=threshold, 
                                 filter=filter_status)

    # plot the trace
    ax.plot(time_axis, trace, color='blue', label=f'trace')

    # add horizontal dashed lines at +/- threshold of the trace
    ax.axhline(y=threshold, color='red', linestyle='--', label=f'+/- threshold={threshold:.2f}')
    ax.axhline(y=-threshold, color='red', linestyle='--')

    # highlight the time windows
    for id, window in enumerate(window_list):
        start_id, stop_id = window
        if id == 0:
            ax.axvspan(start_id, stop_id, color='green', alpha=0.3, label='Time Window(s)')
        else:
            ax.axvspan(start_id, stop_id, color='green', alpha=0.3)

        # annotate time windows
        centre_id = (start_id + stop_id) / 2
        ylim = ax.get_ylim()
        ax.text(centre_id, 1.02*ylim[1], f'[{start_id:.1f}, {stop_id:.1f}]', color='black', fontsize=10, ha='center', va='bottom')

    # set title on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel}, filter {filter_status}', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    # set X-limits
    ax.set_xlim(-4, 4100)

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
    # select 1 good DU to check
    good_du = 1010

    # check 1 ROOT file to test the search method
    check('data/20231028/GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231028001953.141_dat.root', good_du)
    pass

if __name__ == "__main__":
    with record_run_time():
        main()