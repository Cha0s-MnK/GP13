#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

##########################
# PLOT MEAN FFT FOR A DU #
##########################

def load_fft():
    # get information from this file
    basename     = os.path.basename(file)
    du           = int(basename.split('_')[0][2:])
    npz_file     = np.load(file, allow_pickle=True)
    windows_list = {channel: npz_file[f'windows_list{channel}'] for channel in channels}

    # convert GPS times to UTCs
    utcs = gps2utc(npz_file['gps_times'])

    # compute transient rates for 3 ADC channels
    for channel in channels:
        part_hours, part_rates = compute_rate(utcs, windows_list[channel])
        hours_list[channel][du].append(part_hours)
        rates_list[channel][du].append(part_rates)

def plot_fft():
    # create a figure and adjust the figure size
    fig, ax = plt.subplots(figsize=(10,7), dpi=256)

    # plot the FFTs for X (north-south), Y (east-west), Z (up-down)
    colours = ['green', 'blue', 'red']
    for id, channel in enumerate(channels):
        ax.plot(fft_frequency, mean_fft[channel], label=f'Channel {channel}', color=colours[id-1])

    # add vertical lines at 30MHz and 50MHz to test the filter
    ax.axvline(x=30, color='purple', linestyle='--', linewidth=2, label='30 MHz')
    ax.axvline(x=50, color='orange', linestyle='--', linewidth=2, label='50 MHz')

    # set axis scales
    ax.set_yscale('log')

    # set axis labels
    ax.set_xlabel('Frequency / MHz')
    ax.set_ylabel('Mean FFT / a.u.')

    # set figure title
    ax.set_title(f'Mean FFT of DU{du} in {num_entries_du} Entries/Events')

    # set figure legend
    ax.legend(frameon=True)

    # save the figure as a PNG file
    plot_dir  = 'plot/fft/'
    plot_file = f'fft_DU{du}_{num_entries_du}entries_{cutoff_frequency}filter.png'
    plt.savefig(os.path.join(plot_dir, plot_file), bbox_inches='tight')
    print(f'Saved: {plot_dir}{plot_file}')

#################
# MAIN FUNCTION #
#################

def main():
    fft_plot_file = 'result/fft/*.npz'
    
    # get NPZ files and DUs
    file_list, du_list = get_npz_du(fft_plot_file)

    # make dictionaries for easier indexing
    hours_list = {channel: {} for channel in channels}
    rates_list = {channel: {} for channel in channels}

    # loop through all DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            hours_list[channel][du] = []
            rates_list[channel][du] = []

    # loop through DUs to plot FFTs
    for du in du_list:
        plot_fft()

    pass

if __name__ == "__main__":
    with record_run_time():
        main()