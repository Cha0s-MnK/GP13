#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#############################
# PLOT AMPLITUDE HISTOGRAMS #
#############################

def plot_du_amplitude(file_list):
    # get information from the filenames
    du = int(os.path.basename(file_list[0]).split('_')[1][2:])
    date_time, datetime_flat = get_npz_datetime(file_list[0])

    # create a 3x2 layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(3, 2, figsize=(28, 21), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # make dictionaries for easier indexing and initiate them with empty lists
    amplitude_lists = {channel: [] for channel in channels}
    filtered_amplitude_lists = {channel: [] for channel in channels}

    # loop through used files
    for file in file_list:
        # load corresponding NPZ file and get information
        npz_file = np.load(file, allow_pickle=True)
        for channel in channels:
            amplitude_lists[channel].extend(npz_file[f'amplitude_lists{channel}'])
            filtered_amplitude_lists[channel].extend(npz_file[f'filtered_amplitude_lists{channel}'])
    
    # plot a amplitude histogram for each channel
    for id, channel in enumerate(channels):
        plot_channel(axes[2*id], amplitude_lists[channel])
        plot_channel(axes[2*id+1], filtered_amplitude_lists[channel])

        # set channels on the right-hand side
        axes[2*id+1].text(1.01, 0.5, f'channel {channel}', verticalalignment='center', horizontalalignment='left', transform=axes[2*id+1].transAxes, fontsize=18, rotation=-90)

    # set common labels and title
    fig.text(0.5, 0.04, 'Amplitude / ADC units', ha='center', va='center', fontsize=18)
    fig.text(0.04, 0.5, 'Frequency / Hz', ha='center', va='center', rotation='vertical', fontsize=18)
    fig.text(0.25, 0.95, 'No Filter', ha='center', fontsize=18)
    fig.text(0.75, 0.95, f'Filter with Cutoff Frequency = {cutoff_frequency} MHz', ha='center', fontsize=18)
    fig.suptitle(f'Amplitude Distribution for DU{du} on {date_time}'
                 f'\nSample frequency = {sample_frequency} MHz', 
                 fontsize=20)

    # adjust layout: left, bottom, right, top
    plt.tight_layout(rect=[0.04, 0.04, 1.0, 0.98]) 

    # save the figure as a PNG file
    amplitude_plot_file = os.path.join(amplitude_plot_dir, f'amplitude_DU{du}_RUN{num_run}_{datetime_flat}.png')
    plt.savefig(amplitude_plot_file)
    print(f'Saved: {amplitude_plot_file}')

    # close the figure to free up memory
    plt.close()

    pass

def plot_channel(ax,
                 amplitude_list):
    # determine the bins with a width of 1
    width = 1
    bins = np.arange(min(amplitude_list), max(amplitude_list) + 1, width)

    # compute the histogram data
    amplitude_nums, bins, ignored = ax.hist(amplitude_list, bins=bins)
    
    # convert the histogram counts of amplitudes to frequencies
    amplitude_frequencies = amplitude_nums / len(amplitude_list) * 5e8 / width
    
    # re-plot the histogram with frequency values
    ax.clear()
    ax.bar(bins[:-1], amplitude_frequencies, width=width, color='blue', alpha=0.6, label='amplitude distribution')
    
    filtered_amplitude_list = [amp for amp in amplitude_list if -50 <= amp <= 50]
    # compute the Gaussian distribution and convert to frequencies
    amplitude_array = np.array(amplitude_list)
    sorted_amplitude_array = np.sort(amplitude_array)
    low_percentile_value = np.percentile(sorted_amplitude_array, 0.8)
    high_percentile_value = np.percentile(sorted_amplitude_array, 99.2)
    low_id = np.searchsorted(sorted_amplitude_array, low_percentile_value, 'left')
    high_id = np.searchsorted(sorted_amplitude_array, high_percentile_value, 'right')
    Gaussian_amplitude_array = sorted_amplitude_array[low_id:high_id]
    mu = np.mean(Gaussian_amplitude_array)
    sigma = np.std(Gaussian_amplitude_array)
    normal_distribution = norm.pdf(bins, mu, sigma)
    frequency_normal_distribution = normal_distribution * 5e8 / width
    
    # plot the Gaussian distribution
    ax.plot(bins, frequency_normal_distribution, linewidth=2, color='red', label='Gaussian')

    # add vertical lines
    ax.axvline(x=3*sigma, color='orange', linestyle=':', linewidth=2, label='3 sigma')
    ax.axvline(x=-3*sigma, color='orange', linestyle=':', linewidth=2)
    ax.axvline(x=5*sigma, color='green', linestyle=':', linewidth=2, label='5 sigma')
    ax.axvline(x=-5*sigma, color='green', linestyle=':', linewidth=2)

    # set Y-scale
    ax.set_yscale('log')

    # set X-limits and Y-limits
    ax.set_ylim([1e-1, 1e8])

    # enable the grid and legend
    ax.grid(True)
    ax.legend(frameon=True, fontsize=16)

    pass

#################
# MAIN FUNCTION #
#################

def main():
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(amplitude_result_dir, f'amplitude_DU10*_RUN{num_run}*.npz')))

    # get used DUs
    du_list = get_npz_dus(file_list)

    # loop through used DUs
    for du in du_list:
        # get NPZ files for this DU
        file_list = sorted(glob(os.path.join(amplitude_result_dir, f'amplitude_DU{du}_RUN{num_run}*.npz')))

        # plot amplitude histograms for this DU
        plot_du_amplitude(file_list)

    pass

if __name__ == "__main__":
    with record_run_time():
        main()