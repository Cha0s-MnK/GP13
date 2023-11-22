#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#############################
# PLOT MEAN FFTS FOR A FILE #
#############################

def plot_fft2(day_file, 
              night_file):
    # get DU of these files
    basename = os.path.basename(day_file)
    du       = int(basename.split('_')[1][2:])

    # create a 3x2 layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(3, 2, figsize=(36, 24), dpi=200)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the original FFTs and filtered FFTs for each channel in the subplots
    # X(north-south), Y(east-west) and Z(up-down)
    plot_file2(axes, day_file, 'day')
    plot_file2(axes, night_file, 'night')

    # add galactic noise simulations
    galaxy_frequency = range(10, 250) # [MHz]
    for id in range(len(channels)):
        galaxy_noises = np.load(os.path.join(galaxy_sim_dir, galaxy_sim_name[id]))
        # take LST and correct for 20dB linear gain
        galaxy_noise_day   = galaxy_noises[:,14] / linear_gain / linear_gain
        galaxy_noise_night = galaxy_noises[:,2] / linear_gain / linear_gain
        axes[2*id].semilogy(galaxy_frequency, galaxy_noise_day[0:240], "-.", label='Galactic Noise in Day')
        axes[2*id].semilogy(galaxy_frequency, galaxy_noise_night[0:240], "-.", label='Galactic Noise at Night')
        axes[2*id+1].semilogy(galaxy_frequency, galaxy_noise_day[0:240], "-.", label='Galactic Noise in Day')
        axes[2*id+1].semilogy(galaxy_frequency, galaxy_noise_night[0:240], "-.", label='Galactic Noise at Night')

    for id, ax in enumerate(axes):
        # set Y-scale
        ax.set_yscale('log')

        # apply a cutoff by setting the lower limit of the y-axis
        if id % 2 == 0:
            ax.set_ylim(bottom=5e-13)
        else:
            ax.set_ylim(bottom=1e-15)

        # add vertical lines at 30MHz and 50MHz
        ax.axvline(x=30, color='purple', linestyle=':', linewidth=2, label='30 MHz')
        ax.axvline(x=50, color='green', linestyle=':', linewidth=2, label='50 MHz')

        # hide X-ticks of above subplots
        if id != 4 and id != 5:
            ax.set_xticks([])

        # add legends
        ax.legend(frameon=True)

    # set common labels and title
    fig.text(0.5, 0.0, 'Frequency / MHz', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'Mean FFT / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
    fig.text(0.25, 0.95, 'No Filter', ha='center', fontsize=18)
    fig.text(0.75, 0.95, f'Filter with Cutoff Frequency = {cutoff_frequency} MHz', ha='center', fontsize=18)
    fig.suptitle(f'Mean FFT Analysis for DU{du}'
                 f'\nSample frequency = {sample_frequency} MHz', 
                 fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 1, 0.98])

    # save the figure as a PNG file
    fft_plot_name = f'fft_DU{du}_frequency{sample_frequency}_cutoff{cutoff_frequency}.png'
    fft_plot_file = os.path.join(fft_plot_dir, fft_plot_name)
    plt.savefig(fft_plot_file)
    print(f'Saved: {fft_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

    pass

def plot_file2(axes,
              file,
              time_status):
    # get information of this file
    basename                 = os.path.basename(file)
    date_time, datetime_flat = get_npz_datetime(basename)
    npz_file                 = np.load(file, allow_pickle=True)
    mean_fft                 = {channel: npz_file[f'mean_fft{channel}'] for channel in channels}
    mean_filtered_fft        = {channel: npz_file[f'mean_filtered_fft{channel}'] for channel in channels}

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
        axes[2*id].plot(fft_frequency, mean_fft[channel], label=f'{time_status}: {date_time}', color=colour)

        # filtered traces will be on the last three subplots: 1, 3, 5
        axes[2*id+1].plot(fft_frequency, mean_filtered_fft[channel], label=f'{time_status}: {date_time}', color=colour)

        # set the title on the right-hand side
        axes[2*id+1].text(1.01, 0.5, f'channel {channel}', verticalalignment='center', horizontalalignment='left', transform=axes[2*id+1].transAxes, fontsize=16, rotation=-90)
    
    pass

#################
# MAIN FUNCTION #
#################

# plot FFTs of DU1076
def main():
    # get NPZ files containing FFTs
    day_file   = 'result/fft/DU1076/fft_DU1076_frequency500_cutoff50_20231120140803.npz'
    night_file = 'result/fft/DU1076/fft_DU1076_frequency500_cutoff50_20231121013803.npz'

    # plot FFTs and compare them with galactic noises
    plot_fft2(day_file, night_file)

    pass

if __name__ == "__main__":
    with record_run_time():
        main()