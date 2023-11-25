#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#############################
# PLOT MEAN FFTS FOR A FILE #
#############################

def plot_fft(before_day_file, 
             before_night_file, 
             after_day_file, 
             after_night_file):
    # get DU of these files
    basename = os.path.basename(before_day_file)
    du       = int(basename.split('_')[1][2:])

    # create a 3x2 layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(3, 2, figsize=(36, 24), dpi=200)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the original FFTs and filtered FFTs for each channel in the subplots
    # X(north-south), Y(east-west) and Z(up-down)
    plot_file(axes, before_day_file, 'before', 'day')
    plot_file(axes, before_night_file, 'before', 'night')
    plot_file(axes, after_day_file, 'after', 'day')
    plot_file(axes, after_night_file, 'after', 'night')

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
        ax.set_ylim(bottom=1e-12)

        # add vertical lines at 30MHz and 50MHz
        ax.axvline(x=30, color='purple', linestyle=':', linewidth=2, label='30 MHz')
        ax.axvline(x=50, color='orange', linestyle=':', linewidth=2, label='50 MHz')

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
        axes[2*id].plot(fft_frequency, mean_fft[channel], label=f'{time_status} {date_status} trip: {date_time}', linestyle=linestyle, color=colour)

        # filtered traces will be on the last three subplots: 1, 3, 5
        axes[2*id+1].plot(fft_frequency, mean_filtered_fft[channel], label=f'{time_status} {date_status} trip: {date_time}', linestyle=linestyle, color=colour)

        # set the title on the right-hand side
        axes[2*id+1].text(1.01, 0.5, f'channel {channel}', verticalalignment='center', horizontalalignment='left', transform=axes[2*id+1].transAxes, fontsize=16, rotation=-90)
    
    pass

#################
# MAIN FUNCTION #
#################

def main(): 
    # loop through each good DU
    for du in good_du_list:
        # construct filenames for each DU
        before_day_file   = os.path.join(fft_result_dir, f'DU{du}', f'fft_DU{du}_frequency500_cutoff50_20231015141013.npz')
        before_night_file = os.path.join(fft_result_dir, f'DU{du}', f'fft_DU{du}_frequency500_cutoff50_20231014015122.npz')
        after_day_file    = os.path.join(fft_result_dir, f'DU{du}', f'fft_DU{du}_frequency500_cutoff50_20231028140424.npz')
        after_night_file  = os.path.join(fft_result_dir, f'DU{du}', f'fft_DU{du}_frequency500_cutoff50_20231029020123.npz')

        # plot and compare FFTs before and after the site trip for each DU
        plot_fft(before_day_file, before_night_file, after_day_file, after_night_file)

    pass

if __name__ == "__main__":
    with record_run_time():
        main()