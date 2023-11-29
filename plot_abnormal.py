#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

###################################
# PLOT ABNORMAL TRACE AND ITS PSD #
###################################

def plot_abnormal(file):
    # get used DU of this file from the filename
    du = int(os.path.basename(file).split('_')[1][2:])

    # get the date and time of this file
    date_time, datetime_flat = get_npz_datetime(os.path.basename(file))

    # get inner information of this NPZ file
    npz_file   = np.load(file, allow_pickle=True)
    trace      = npz_file['trace']
    channel    = str(npz_file['channel'])
    psd        = npz_file['psd']
    filter_psd = npz_file['filter_psd']

    # filter the trace
    filter_trace = high_pass_filter(trace)

    # create a figure and a GridSpec layout
    fig = plt.figure(figsize=(32, 32), dpi=256)
    gs = gridspec(6, 1)

    # create the subplots
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1:3, 0])
    ax2 = fig.add_subplot(gs[3, 0])
    ax3 = fig.add_subplot(gs[4:, 0])

    # plot the abnormal trace and filtered one
    plot_trace(ax=ax0, du=du, trace=trace, channel=channel)
    plot_trace(ax=ax2, du=du, trace=trace, channel=channel, filter='on')

    # plot corresponding PSDs
    hour = date_time.hour
    plot_psd(ax=ax1, channel=channel, hour=hour, psd=psd)
    plot_psd(ax=ax3, channel=channel, hour=hour, psd=filter_psd)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 1, 0.98])

    # save the figure as a PNG file
    abnormal_plot_file = os.path.join(abnormal_plot_dir, f'abnormal_DU{du}_RUN{num_run}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{datetime_flat}.png')
    plt.savefig(abnormal_plot_file)
    print(f'Saved: {abnormal_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

def plot_trace(ax,
               du,
               trace,
               channel,
               filter='off'):
    # apply the high-pass filter if filter is set to 'on'
    if filter == 'on':
        trace = high_pass_filter(trace)
    
    # plot the trace
    ax.plot(time_axis, trace, color='blue', label=f'trace')

    # set labels and title
    ax.set_xlabel('ADC units', fontsize=16)
    ax.set_ylabel('Time / ns', fontsize=16)
    ax.set_title(f'channel {channel}, filter {filter}', fontsize=18)

    # add horizontal dashed lines at +/- threshold of the trace
    threshold = num_threshold * noises[channel][str(du)]
    ax.axhline(y=threshold, color='red', linestyle='--', label=f'+/- threshold={threshold}')
    ax.axhline(y=-threshold, color='red', linestyle='--')

    # add horizontal dashed lines at +/- STD of the trace
    std = np.std(trace)
    ax.axhline(y=std, color='orange', linestyle='--', label=f'+/- STD={std:.2f}')
    ax.axhline(y=-std, color='orange', linestyle='--')

    # set X-limit
    ax.set_xlim([0, num_samples*time_step])

    # set X-ticks every 250 nanoseconds
    ax.set_xticks(np.arange(min(time_axis), max(time_axis)+1, 250))

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, fontsize=16)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

    pass

def plot_psd(ax,
             channel,
             hour,
             psd):
    # plot the PSD
    ax.plot(fft_frequency, psd, label='PSD', color='red')

    # load galactic noise simulation for this channel
    galaxy_sim = np.load(os.path.join(galaxy_dir, galaxy_name[mapping[channel]]))

    # apply noise according to the hour and correct the 20dB linear gain
    galaxy_noise = galaxy_sim[:,hour] / linear_gain / linear_gain

    # plot the galactic noise simulation
    ax.semilogy(range(10, 250), galaxy_noise[0:240], color='orange', linestyle='-.', label='Galactic Noise')

    # set X-limit
    ax.set_xlim([0, sample_frequency/2])

    # set Y-scale and Y-limit
    ax.set_yscale('log')
    #ax.set_ylim(bottom=1e-13)

    # add vertical lines at 30MHz and 50MHz
    ax.axvline(x=30, color='purple', linestyle=':', linewidth=2, label='30 MHz')
    ax.axvline(x=50, color='green', linestyle=':', linewidth=2, label='50 MHz')

    # set X-ticks every 25MHz
    ax.xaxis.set_major_locator(ticker.MultipleLocator(25))

    # set labels and title
    ax.set_xlabel('Frequency / MHz', fontsize=16)
    ax.set_ylabel('Mean FFT PSD / $V^2 MHz^{-1}$', fontsize=16)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, fontsize=16)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

    pass

#################
# MAIN FUNCTION #
#################

# plot mean FFT PSDs for each DU
def main():
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(abnormal_result_dir, '*.npz')))

    # loop through all files
    for file in file_list:
        plot_abnormal(file)

    pass

if __name__ == "__main__":
    with record_run_time():
        main()