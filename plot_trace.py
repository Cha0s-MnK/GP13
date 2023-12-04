#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

# create the save directory if it does not exist
save_dir = os.path.join(plot_dir, f'{trace_wanted}-trace')
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

# plot the traces and corresponding PSDs for 1 entry
def plot1npz_trace(file):
    # get information of this file
    du                  = int(os.path.basename(file).split('_')[2][2:])
    trace_wanted        = os.path.basename(file).split('-')[0]
    DHtime, DHtime_flat = get_npz_DHtime(os.path.basename(file))
    hour                = DHtime.hour
    npz_file            = np.load(file, allow_pickle=True)
    traces              = {channel: npz_file[f'traces{channel}'] for channel in channels}
    filter_traces       = {channel: npz_file[f'filter_traces{channel}'] for channel in channels}
    psds                = {channel: npz_file[f'psds{channel}'] for channel in channels}
    filter_psds         = {channel: npz_file[f'filter_psds{channel}'] for channel in channels}

    # adjust the figure size and create a GridSpec layout for the subplots
    fig = plt.figure(figsize=(60, 32), dpi=165)
    gs  = gridspec(12, 2)

    # make dictionaries for easier indexing
    axes   = {}
    layout = {
        '01': (0, 0, 1, 1),
        '02': (1, 0, 3, 1),
        '03': (4, 0, 1, 1),
        '04': (5, 0, 3, 1),
        '05': (8, 0, 1, 1),
        '06': (9, 0, 3, 1),
        '07': (0, 1, 1, 1),
        '08': (1, 1, 3, 1),
        '09': (4, 1, 1, 1),
        '10': (5, 1, 3, 1),
        '11': (8, 1, 1, 1),
        '12': (9, 1, 3, 1),
    }

    # loop through the layout dictionary to create each subplot
    for ax_id, (row, column, width, height) in layout.items():
        axes[ax_id] = fig.add_subplot(gs[row:row+width, column:column+height])

    # loop through 3 ADC channels to plot the traces, the filtered ones and corresponding PSDs
    for i, channel in enumerate(channels):
        base_id = i * 2
        subplot_trace(ax=axes[f'{base_id + 1:02}'], channel=channel, trace=traces[channel], du=du, filter_status='off')
        subplot_psd(ax=axes[f'{base_id + 2:02}'], channel=channel, psd=psds[channel], hour=hour)
        subplot_trace(ax=axes[f'{base_id + 7:02}'], channel=channel, trace=filter_traces[channel], du=du, filter_status='on')
        subplot_psd(ax=axes[f'{base_id + 8:02}'], channel=channel, psd=filter_psds[channel], hour=hour)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 1, 0.96])

    # set common labels and title
    #fig.text(0.5, 0.0, 'Frequency / MHz', ha='center', fontsize=18)
    #fig.text(0.0, 0.5, 'Mean FFT / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
    fig.suptitle(f'Traces, Filtered Traces and Corresponding PSDs\nRUN{num_run}, DU{du}, Time = {DHtime}, Sample frequency = {sample_frequency} MHz, Cut-off Frequency = {cutoff_frequency} MHz', fontsize=20)

    # save the figure as a PNG file
    save_file = os.path.join(save_dir, f'{trace_wanted}-trace_RUN{num_run}_DU{du}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{DHtime_flat}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')

    # close the figure to free up memory
    plt.close(fig)

# plot 1 trace on 1 ax
def subplot_trace(ax, channel, trace, du, filter_status):
    # plot the trace
    ax.plot(time_axis, trace, color='blue', label=f'trace')

    # add horizontal dashed lines at +/-STD of the trace
    std = np.std(trace)
    ax.axhline(y=std, color='orange', linestyle='--', label=f'+/- STD={std:.2f}')
    ax.axhline(y=-std, color='orange', linestyle='--')

    # add horizontal dashed lines at +/-threshold of the trace
    threshold = num_threshold * std
    ax.axhline(y=threshold, color='red', linestyle='--', label=f'+/- threshold={threshold:.2f}')
    ax.axhline(y=-threshold, color='red', linestyle='--')

    if trace_wanted == 'pulse' and filter_status == 'on':
        # highlight the transient/pulse windows
        window_list = search_windows_test(trace=trace, num_threshold=num_threshold, filter_status='off')
        for window_id, window in enumerate(window_list):
            start_id, stop_id = window
            start_time        = start_id * time_step
            stop_time         = stop_id * time_step
            if window_id == 0:
                ax.axvspan(start_time, stop_time, color='green', alpha=0.3, label='transient/pulse window(s)')
            else:
                ax.axvspan(start_time, stop_time, color='green', alpha=0.3)

            # annotate time windows
            centre_time = (start_time + stop_time) / 2
            ylim = ax.get_ylim()
            ax.text(centre_time, 1.02*ylim[1], f'[{start_time}, {stop_time}]', color='black', fontsize=14, ha='center', va='bottom')

    # set X-ticks
    ax.set_xlim([min(time_axis), max(time_axis)])
    ax.set_xticks(np.arange(min(time_axis), max(time_axis)+1, 200))

    # set labels and title
    ax.set_xlabel('Time / ns', fontsize=16)
    ax.set_ylabel('ADC units', fontsize=16)
    ax.set_title(f'channel {channel} with filter {filter_status}', fontsize=18, pad=20)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

# plot 1 PSD on 1 ax
def subplot_psd(ax, channel, psd, hour):
    # plot the PSD
    ax.plot(fft_frequency, psd, label='PSD', color='red')

    # load galactic noise simulation for this channel
    galaxy_sim = np.load(os.path.join(galaxy_dir, galaxy_name[channel_mapping[channel]]))

    # apply noise according to the hour and correct the 20dB linear gain
    galaxy_noise = galaxy_sim[:,hour] / linear_gain / linear_gain

    # plot the galactic noise simulation
    ax.semilogy(range(10, 250), galaxy_noise[0:240], color='orange', linestyle='-.', label='Galactic Noise')

    # add vertical lines at 30MHz and 50MHz
    ax.axvline(x=30, color='purple', linestyle=':', linewidth=2, label='30 MHz')
    ax.axvline(x=50, color='green', linestyle=':', linewidth=2, label='50 MHz')

    # set X-ticks
    ax.set_xlim([min(fft_frequency), max(fft_frequency)])
    ax.set_xticks(np.arange(min(fft_frequency), max(fft_frequency)+1, 25))

    # set Y-scale and Y-limit
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-14)

    # set labels and title
    ax.set_xlabel('Frequency / MHz', fontsize=16)
    ax.set_ylabel('PSD / $V^2 MHz^{-1}$', fontsize=16)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

def get_mean_psd(result_dir):
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(result_dir, '*.npz')))

    # get used DUs
    du_list = get_npz_dus(file_list)

    # loop through DUs
    for du in du_list:
        mean_psd, mean_filter_psd, num_files = 0, 0, 0
        # loop through files
        for file in file_list:
            # check used DU of this file
            if du != int(os.path.basename(file).split('_')[1][2:]):
                continue

            # check channel of this file
            npz_file   = np.load(file, allow_pickle=True)
            trace      = npz_file['trace']
            channel    = str(npz_file['channel'])
            if channel != 'X':
                continue

            # get information of this NPZ file
            trace      = npz_file['trace']

            # get FFT PSD and filtered FFT PSD
            psd        = get_psd(traces[channel_id])
            filter_psd = get_psd(filter_trace)

            # add the PSD
            mean_psd        += psd
            mean_filter_psd += filter_psd
            num_files       += 1

        # compute the mean PSD
        mean_psd        /= num_files
        mean_filter_psd /= num_files

        # create a layout for the subplots and adjust the figure size
        fig, axes = plt.subplots(2, 1, figsize=(20, 16), dpi=256)

        # flatten axes array for easier indexing
        axes = axes.flatten()

        # plot the mean PSD and the filtered one
        channel = 'X'
        hour = 8 # mean hour?
        plot_psd(ax=axes[0], psd=mean_psd, channel=channel, hour=hour)
        plot_psd(ax=axes[1], psd=mean_filter_psd, channel=channel, hour=hour)

        # set common labels and title
        #fig.text(0.5, 0.0, 'Duration Time / ns', ha='center', fontsize=18)
        #fig.text(0.0, 0.5, 'Probability*20', va='center', rotation='vertical', fontsize=18)
        plt.suptitle(f'Mean FFT PSD for DU{du} channel {channel} in RUN{num_run}', fontsize=20)

        # adjust the layout: left, bottom, right, top
        plt.tight_layout(rect=[0.01, 0.01, 1, 0.98])

        # save the figure as a PNG file
        abnormal_plot_file = os.path.join(abnormal_plot_dir, f'abnormal-mean-PSD_DU{du}_RUN{num_run}_channel{channel}_frequency{sample_frequency}_cutoff{cutoff_frequency}.png')
        plt.savefig(abnormal_plot_file)
        print(f'Saved: {abnormal_plot_file}')

        # close the figure to free up memory
        plt.close(fig)

#################
# MAIN FUNCTION #
#################

# plot mean FFT PSDs for each DU
def main():
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(result_dir, f'{trace_wanted}-trace', f'{trace_wanted}-trace*RUN{num_run}*.npz')))
    num_files = len(file_list)

    # get used DUs
    du_list = get_npz_dus(file_list=file_list)

    # loop through listed files
    print(f'\nLoop through {num_files} NPZ files from RUN{num_run}.\n')
    for file in file_list:
        plot1npz_trace(file)

if __name__ == "__main__":
    with record_run_time():
        main()