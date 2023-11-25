#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

##############################
# PLOT MEAN FFT PSD FOR A DU #
##############################

def plot_fft(day_file, 
             night_file):
    # get used DU of these files from the filename
    du = int(os.path.basename(day_file).split('_')[1][2:])

    # get the day and night time
    day, day_flat     = get_npz_datetime(os.path.basename(day_file))
    night, night_flat = get_npz_datetime(os.path.basename(day_file))

    # create a 3x2 layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(3, 2, figsize=(36, 24), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the original FFTs and filtered FFTs for each channel in the subplots
    plot_file(axes, day_file, 'day')
    plot_file(axes, night_file, 'night')

    # add galactic noise simulations
    galaxy_frequency = range(10, 250) # [MHz]
    for id, channel in enumerate(channels):
        # load noise simulation for this channel
        galaxy_noises = np.load(os.path.join(galaxy_dir, galaxy_name[id]))

        # take noise according to the time and correct for 20dB linear gain
        galaxy_noise_day   = galaxy_noises[:,14] / linear_gain / linear_gain
        galaxy_noise_night = galaxy_noises[:,2] / linear_gain / linear_gain

        # plot galactic noise simulations
        axes[2*id].semilogy(galaxy_frequency, galaxy_noise_day[0:240], color='orange', linestyle='-.', label='Day Galactic Noise')
        axes[2*id].semilogy(galaxy_frequency, galaxy_noise_night[0:240], color='cyan', linestyle='-.', label='Night Galactic Noise')
        axes[2*id+1].semilogy(galaxy_frequency, galaxy_noise_day[0:240], color='orange', linestyle='-.', label='Day Galactic Noise')
        axes[2*id+1].semilogy(galaxy_frequency, galaxy_noise_night[0:240], color='cyan', linestyle='-.', label='Night Galactic Noise')

        # set comments on the right-hand side
        axes[2*id+1].text(1.01, 0.5, f'channel {channel}', verticalalignment='center', horizontalalignment='left', transform=axes[2*id+1].transAxes, fontsize=18, rotation=-90)

    for ax in axes:
        # set Y-scale
        ax.set_yscale('log')

        # apply a cutoff by setting the lower limit of the y-axis
        ax.set_ylim(bottom=1e-13)

        # add vertical lines at 30MHz and 50MHz
        ax.axvline(x=30, color='purple', linestyle=':', linewidth=2, label='30 MHz')
        ax.axvline(x=50, color='green', linestyle=':', linewidth=2, label='50 MHz')

        # set X-ticks every 25MHz
        ax.xaxis.set_major_locator(ticker.MultipleLocator(25))

        # enable the grid and legend
        ax.grid(True)
        ax.legend(frameon=True, fontsize=16)

    # set common labels and title
    fig.text(0.5, 0.0, 'Frequency / MHz', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'Mean FFT / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
    fig.text(0.25, 0.95, 'No Filter', ha='center', fontsize=18)
    fig.text(0.75, 0.95, f'Filter with Cutoff Frequency = {cutoff_frequency} MHz', ha='center', fontsize=18)
    fig.suptitle(f'Mean FFT PSD Analysis for DU{du}'
                 f'\nSample frequency = {sample_frequency} MHz, Day Time = {day}, Night Time = {night}', 
                 fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 1, 0.98])

    # save the figure as a PNG file
    fft_plot_file = os.path.join(fft_plot_dir, f'fft_DU{du}_frequency{sample_frequency}_cutoff{cutoff_frequency}_day{day_flat}_night{night_flat}.png')
    plt.savefig(fft_plot_file)
    print(f'Saved: {fft_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

    pass

def plot_file(axes,
              file,
              time_status):
    # get information of this file
    date_time, datetime_flat = get_npz_datetime(os.path.basename(file))
    npz_file                 = np.load(file, allow_pickle=True)
    mean_fft                 = {channel: npz_file[f'mean_fft{channel}'] for channel in channels}
    mean_filtered_fft        = {channel: npz_file[f'mean_filtered_fft{channel}'] for channel in channels}

    # change colour based on the time status
    if time_status == 'day':
        colour = 'red'
    else:
        colour = 'blue'

    # plot original FFTs and filtered FFTs for each channel in the subplots
    for id, channel in enumerate(channels):
        # original traces will be on left three subplots: 0, 2, 4
        axes[2*id].plot(fft_frequency, mean_fft[channel], label=f'{time_status}', color=colour)

        # filtered traces will be on the last three subplots: 1, 3, 5
        axes[2*id+1].plot(fft_frequency, mean_filtered_fft[channel], label=f'{time_status}', color=colour)
    
    pass

#################
# MAIN FUNCTION #
#################

# plot mean FFT PSDs for each DU
def main():
    # get the string of NPZ files
    files_str = fft_result_dir + '/*/*.npz'

    # get used DUs
    du_list = get_npz_du(files_str)

    # loop through used DUs
    for du in du_list:
        # get NPZ files for each DU
        day_file   = os.path.join(fft_result_dir, f'DU{du}', f'fft_DU{du}_frequency500_cutoff50_20231116140503.npz')
        night_file = os.path.join(fft_result_dir, f'DU{du}', f'fft_DU{du}_frequency500_cutoff50_20231116015403.npz')

        # plot mean FFT PSDs and compare them with galactic noise simulations
        plot_fft(day_file, night_file)

    pass

if __name__ == "__main__":
    with record_run_time():
        main()