#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#############################
# PLOT MEAN FFTS FOR A FILE #
#############################

def plot1fft(file):
    # number of samples in a trace
    num_samples = 1024

    # frequencies of the FFT [MHz]
    fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency

    # get information of this file
    basename        = os.path.basename(file)
    du              = int(basename.split('_')[1][2:])
    date_time, datetime_flat = get_npz_datetime(basename)
    npz_file        = np.load(file, allow_pickle=True)
    mean_fft        = {channel: npz_file[f'mean_fft{channel}'] for channel in channels}
    mean_filter_fft = {channel: npz_file[f'mean_filtered_fft{channel}'] for channel in channels}

    # create a figure with 2 subplots and adjust the figure size
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 14), dpi=256)

    # plot the original FFTs and filtered FFTs for each channel in the subplots
    # X(north-south), Y(east-west) and Z(up-down)
    colours = ['green', 'blue', 'red']
    for id, channel in enumerate(channels):
        ax1.plot(fft_frequency, mean_fft[channel], label=f'Channel {channel}', color=colours[id])
        ax1.set_ylabel('Mean FFT / a.u.')

        ax2.plot(fft_frequency, mean_filter_fft[channel], label=f'Channel {channel}', color=colours[id])
        ax2.set_xlabel('Frequency / MHz')
        ax2.set_ylabel('Mean Filtered FFT / a.u.')

    # set Y-scale
    ax1.set_yscale('log')
    ax2.set_yscale('log')

    # add vertical lines at 30MHz and 50MHz to test the filter in both subplots
    for ax in [ax1, ax2]:
        ax.axvline(x=30, color='purple', linestyle='--', linewidth=2, label='30 MHz')
        ax.axvline(x=50, color='orange', linestyle='--', linewidth=2, label='50 MHz')
    
    # add legends
    ax1.legend(frameon=True)
    ax2.legend(frameon=True)

    # hide the X-ticks of ax1
    ax1.set_xticks([])

    # add a supertitle for the entire figure
    fig.suptitle(f'Mean FFT Analysis for DU{du} on {date_time}'
                 f'\nThreshold = {num_threshold}, Separation = {standard_separation} samples, Least Crossing = {num_crossings}'
                 f'\nMax samples = {max_samples}, Sample frequency = {sample_frequency} MHz, Cutoff frequency = {cutoff_frequency} MHz', 
                 fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    # save the figure as a PNG file
    fft_plot_name = f'fft_DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{datetime_flat}.png'
    fft_plot_file = os.path.join(fft_plot_dir, fft_plot_name)
    plt.savefig(fft_plot_file)
    print(f'Saved: {fft_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

#################
# MAIN FUNCTION #
#################

def main(): 
    # get NPZ files
    file_list = sorted(glob.glob('result/fft_old/*.npz'))

    # loop through NPZ files to plot FFTs
    for file in file_list:
        plot1fft(file)

    pass

if __name__ == "__main__":
    with record_run_time():
        main()