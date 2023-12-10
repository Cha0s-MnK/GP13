#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

# create the save directory if it does not exist
name = f'{wanted}-trace'
save_dir = os.path.join(plot_dir, name)
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

# plot the traces and corresponding PSDs for 1 entry
def plot1npz_trace(file):
    # get info from this file
    du            = get_npz_du(file=file)
    cst, cst_flat = get_npz_cst(file=file)
    hour          = cst.hour
    npz_file      = np.load(file, allow_pickle=True)
    traces        = {channel: npz_file[f'traces{channel}'] for channel in channels}
    filter_traces = {channel: npz_file[f'filter_traces{channel}'] for channel in channels}
    psds          = {channel: npz_file[f'psds{channel}'] for channel in channels}
    filter_psds   = {channel: npz_file[f'filter_psds{channel}'] for channel in channels}

    # adjust the figure size and create a GridSpec layout for the subplots
    fig = plt.figure(figsize=(32, 45), dpi=165)
    gs  = gridspec(12, 2)

    # make dictionaries for easier indexing
    axes   = {}
    layout = {
        '01': (0, 0, 1, 2),
        '02': (2, 0, 2, 1),
        '03': (1, 0, 1, 2),
        '04': (2, 1, 2, 1),
        '05': (4, 0, 1, 2),
        '06': (6, 0, 2, 1),
        '07': (5, 0, 1, 2),
        '08': (6, 1, 2, 1),
        '09': (8, 0, 1, 2),
        '10': (10, 0, 2, 1),
        '11': (9, 0, 1, 2),
        '12': (10, 1, 2, 1),
    }

    # loop through the layout dictionary to create each subplot
    for ax_id, (row, column, width, height) in layout.items():
        axes[ax_id] = fig.add_subplot(gs[row:row+width, column:column+height])

    # loop through 3 ADC channels to plot the traces, the filtered ones and corresponding PSDs
    for i, channel in enumerate(channels):
        base_id = i * 4
        subplot_trace(ax=axes[f'{base_id + 1:02}'], channel=channel, trace=traces[channel], du=du, filter_status='off')
        subplot_psd(ax=axes[f'{base_id + 2:02}'], channel=channel, psd=psds[channel], hour=hour, filter_status='off')
        subplot_trace(ax=axes[f'{base_id + 3:02}'], channel=channel, trace=filter_traces[channel], du=du, filter_status='on')
        subplot_psd(ax=axes[f'{base_id + 4:02}'], channel=channel, psd=filter_psds[channel], hour=hour, filter_status='on')

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 1, 0.96])

    # set common labels and title
    #fig.text(0.5, 0.0, 'Frequency / MHz', ha='center', fontsize=18)
    #fig.text(0.0, 0.5, 'Mean FFT / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
    fig.suptitle(f'Traces, Filtered Traces and Corresponding PSDs\nRUN{nums_run}, DU{du}, Time = {cst}, Cut-off Frequency = {cutoff_frequency} MHz', fontsize=20)

    # save the figure as a PNG file
    save_file = os.path.join(save_dir, f'{name}_RUN{nums_run}_DU{du}_cutoff{cutoff_frequency}_{cst_flat}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')

    plt.close(fig)

# plot 1 trace on 1 ax
def subplot_trace(ax, channel, trace, du, filter_status):
    # plot the trace
    ax.plot(time_axis, trace, color='blue')

    # add horizontal dashed lines at +/-thresholds of the trace
    #threshold1 = num_threshold1*rms(trace)
    #threshold2 = num_threshold2*rms(trace)
    ax.axhline(y=threshold1, color='red', linestyle='--', label=f'+/-Threshold1={threshold1:.2f}')
    ax.axhline(y=-threshold1, color='red', linestyle='--')
    ax.axhline(y=threshold2, color='darkorange', linestyle='--', label=f'+/-Threshold2={threshold2:.2f}')
    ax.axhline(y=-threshold2, color='darkorange', linestyle='--')

    # highlight the transient/pulse windows
    window_list = search_windows(trace=trace, filter_status='off')
    for window_id, window in enumerate(window_list):
        start_id, stop_id = window
        start_time        = start_id * time_step
        stop_time         = stop_id * time_step
        ax.axvspan(start_time, stop_time, color='green', alpha=0.3)

        # annotate time windows
        centre_time = (start_time + stop_time) / 2
        ylim = ax.get_ylim()
        ax.text(centre_time, 1.02*ylim[1], f'[{start_time}, {stop_time}]', color='black', fontsize=14, ha='center', va='bottom')

    # set X-ticks
    ax.set_xlim([min(time_axis), max(time_axis)])
    ax.set_xticks(np.arange(min(time_axis), max(time_axis)+1, 200))

    # settings
    ax.set_xlabel('Time / ns', fontsize=16)
    ax.set_ylabel('ADC units', fontsize=16)
    ax.set_title(f'channel {channel} with filter {filter_status}', fontsize=18, pad=20)
    ax.grid(True)
    ax.legend(ncol=2, frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

# plot 1 PSD on 1 ax
def subplot_psd(ax, channel, psd, hour, filter_status):
    # plot the PSD
    ax.plot(fft_frequency, psd, label='PSD', color='red')

    # load galactic noise simulation for this channel
    galaxy_sim = np.load(os.path.join(galaxy_dir, galaxy_name[channel_mapping[channel]]))

    # apply noise according to the hour and correct the 20dB linear gain
    galaxy_noise = galaxy_sim[:,hour] / linear_gain / linear_gain

    # plot the galactic noise simulation
    ax.semilogy(range(10, 250), galaxy_noise[0:240], color='orange', linestyle='-.', label='Galactic Noise')

    # add vertical lines at 30MHz and 50MHz
    #ax.axvline(x=30, color='purple', linestyle=':', linewidth=2, label='30 MHz')
    ax.axvline(x=50, color='green', linestyle=':', linewidth=2, label='50 MHz')

    # set X-ticks
    ax.set_xlim([min(fft_frequency), max(fft_frequency)])
    ax.set_xticks(np.arange(min(fft_frequency), max(fft_frequency)+1, 25))

    # set Y-scale and Y-limit
    ax.set_yscale('log')
    ax.set_ylim([1e-14, 1e-6])

    # settings
    ax.set_xlabel('Frequency / MHz', fontsize=16)
    ax.set_ylabel('PSD / $V^2 MHz^{-1}$', fontsize=16)
    ax.set_title(f'channel {channel} with filter {filter_status}', fontsize=18, pad=10)
    ax.grid(True)
    ax.legend(ncol=1, frameon=True, loc='upper right', fontsize=14)
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

def main():
    file_list = get_npz_files(name=name)

    for file in file_list:
        plot1npz_trace(file)

if __name__ == "__main__":
    with record_run_time():
        main()