#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#############################
# PLOT MEAN FFTS FOR A FILE #
#############################

def compare_fft(before_day_file, after_day_file, before_night_file, after_night_file):
    # get DU of these files
    basename = os.path.basename(before_day_file)
    du       = int(basename.split('_')[1][2:])

    # create a 3x2 layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(3, 2, figsize=(32, 20), dpi=200)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the original FFTs and filtered FFTs for each channel in the subplots
    # X(north-south), Y(east-west) and Z(up-down)
    before_day_datetime   = plot_file(axes, before_day_file, 'before', 'day')
    before_night_datetime = plot_file(axes, before_night_file, 'before', 'night')
    after_day_datetime    = plot_file(axes, after_day_file, 'after', 'day')
    after_night_datetime  = plot_file(axes, after_night_file, 'after', 'night')

    # add galactic noise simulations
    galaxy_frequency = range(10, 250) # [MHz]
    for id in range(len(channels)):
        galaxy_noises = np.load(os.path.join(galaxy_sim_dir, galaxy_sim_name[id]))
        galaxy_noise_day   = galaxy_noises[:,13] / linear_gain / linear_gain  #Take 18h LST - correct for 20dB gain
        galaxy_noise_night = galaxy_noises[:,0] / linear_gain / linear_gain   #Take 6h LST - correct for 20dB gain
        axes[2*id].semilogy(galaxy_frequency, galaxy_noise_day[0:240], "--", label='Galaxy Noise in Day (Sandra)')
        axes[2*id].semilogy(galaxy_frequency, galaxy_noise_night[0:240], "--", label='Galaxy Noise at Night (Sandra)')
        axes[2*id+1].semilogy(galaxy_frequency, galaxy_noise_day[0:240], "--", label='Galaxy Noise in Day (Sandra)')
        axes[2*id+1].semilogy(galaxy_frequency, galaxy_noise_night[0:240], "--", label='Galaxy Noise at Night (Sandra)')

    for id, ax in enumerate(axes):
        # set Y-scale
        ax.set_yscale('log')

        # add vertical lines at 30MHz and 50MHz
        ax.axvline(x=30, color='purple', linestyle='--', linewidth=2, label='30 MHz')
        ax.axvline(x=50, color='orange', linestyle='--', linewidth=2, label='50 MHz')

        # hide X-ticks of above subplots
        if id != 4 and id != 5:
            ax.set_xticks([])

        # add legends
        ax.legend(frameon=True)

    # set common labels and title
    fig.text(0.5, 0.0, 'Frequency / MHz', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'Mean FFT / ADC units', va='center', rotation='vertical', fontsize=18)
    fig.suptitle(f'Mean FFT Analysis for DU{du}'
                 f'\nSample frequency = {sample_frequency} MHz, Cutoff frequency = {cutoff_frequency} MHz', 
                 fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 1, 0.98])

    # save the figure as a PNG file
    compare_plot_name = f'compare_DU{du}_frequency{sample_frequency}_cutoff{cutoff_frequency}.png'
    compare_plot_file = os.path.join(compare_plot_dir, compare_plot_name)
    plt.savefig(compare_plot_file)
    print(f'Saved: {compare_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

    pass

def plot_file(axes,
              file,
              date_status,
              time_status):
    # get information of this file
    basename                 = os.path.basename(file)
    date_time, datetime_flat = get_npz_datetime(basename)
    npz_file                 = np.load(file, allow_pickle=True)
    mean_fft                 = {channel: npz_file[f'mean_fft{channel}'] for channel in channels}
    mean_filtered_fft        = {channel: npz_file[f'mean_filtered_fft{channel}'] for channel in channels}

    # get FFT frequency based on date status
    if date_status == 'before':
        # number of samples in a trace
        num_samples = 1024

        # frequencies of the FFT [MHz]
        fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency
    else:
        # number of samples in a trace
        num_samples = 2048

        # frequencies of the FFT [MHz]
        fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency

    # change line style based on date status
    if date_status == 'before':
        linestyle = '--'
    else:
        linestyle = '-'

    # change colour based on time status
    if time_status == 'day':
        colour = 'red'
    else:
        colour = 'blue'

    # plot original FFTs and filtered FFTs for each channel in the subplots
    # X(north-south), Y(east-west) and Z(up-down)
    colours = ['green', 'blue', 'red']
    for id, channel in enumerate(channels):
        # original traces will be on left three subplots: 0, 2, 4
        axes[2*id].plot(fft_frequency, mean_fft[channel], label=f'{date_status} {time_status}', linestyle=linestyle, color=colour)

        # filtered traces will be on the last three subplots: 1, 3, 5
        axes[2*id+1].plot(fft_frequency, mean_filtered_fft[channel], label=f'{date_status} {time_status}', linestyle=linestyle, color=colour)

        # set the title on the right-hand side
        axes[2*id+1].text(1.01, 0.5, f'channel {channel}', verticalalignment='center', horizontalalignment='left', transform=axes[2*id+1].transAxes, fontsize=16, rotation=-90)
    
    return date_time

#################
# MAIN FUNCTION #
#################

def main(): 
    # get NPZ files
    before_day_file   = 'result/fft/DU1010/fft_DU1010_frequency500_cutoff50_20231014131923.npz'
    after_day_file    = 'result/fft/DU1010/fft_DU1010_frequency500_cutoff50_20231028121653.npz'
    before_night_file = 'result/fft/DU1010/fft_DU1010_frequency500_cutoff50_20231014015122.npz'
    after_night_file  = 'result/fft/DU1010/fft_DU1010_frequency500_cutoff50_20231028001953.npz'

    # compare FFTs before and after the site trip
    compare_fft(before_day_file, after_day_file, before_night_file, after_night_file)

    pass

if __name__ == "__main__":
    with record_run_time():
        main()