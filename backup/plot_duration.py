#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

################################################
# COMPUTE AND PLOT TRANSIENT RATES FOR EACH DU #
################################################

def get_duration(du):
    # update NPZ files according to current DU
    file_list = sorted(glob(os.path.join(duration_result_dir, f'*/*DU{du}_RUN{num_run}*.npz')))

    # make dictionaries for easier indexing and initiate them with empty lists
    all_durations_list = {channel: [] for channel in channels}

    # loop through NPZ files
    for file in file_list:
        # get information from this file
        npz_file       = np.load(file, allow_pickle=True)
        durations_list = {channel: npz_file[f'durations_list{channel}'] for channel in channels}

        # gather durations for each channel
        for channel in channels:
            all_durations_list[channel].extend(durations_list[channel])

    # plot the duration time distribution of transients / pulses
    plot_duration(du, all_durations_list)

    pass
    
def plot_duration(du,
                  durations_list):
    # create a layout for the subplots and adjust the figure size
    fig, axes = plt.subplots(3, 1, figsize=(30, 15), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the duration time distribution of transients / pulses for each channel
    for id, channel in enumerate(channels):
        plot_duration_channel(axes[id], channel, durations_list[channel])

    # set common labels and title
    fig.text(0.5, 0.0, 'Duration Time / ns', ha='center', fontsize=18)
    fig.text(0.0, 0.5, 'Probability*20', va='center', rotation='vertical', fontsize=18)
    plt.suptitle(f'Duration Time Distribution of Transients/Pulses for DU{du} in RUN{num_run}', fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])

    # create the save directory if it does not exist
    os.makedirs(duration_plot_dir, exist_ok=True)

    # save the figure as a PNG file
    duration_plot_file = os.path.join(duration_plot_dir, f'duration_DU{du}_RUN{num_run}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_fluctuation{std_fluctuation}_frequency{sample_frequency}_cutoff{cutoff_frequency}.png')
    plt.savefig(duration_plot_file)
    print(f'Saved: {duration_plot_file}')

    # close the figure to free up memory
    plt.close(fig)

def plot_duration_channel(ax, 
                          channel,
                          duration_list):
    if len(duration_list) == 0:
        return

    # determine the bins with a width of 1
    width = 20
    bins = np.arange(0, num_samples*time_step, width)

    # plot the histogram
    ax.hist(duration_list, bins=bins, density=True, color='blue', alpha=0.75, label='duration time distribution')

    # set X-limits and Y-limits
    ax.set_xlim([0, num_samples*time_step])
    #ax.set_ylim(bottom=1e-1)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', fontsize=16)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

    # set the channel on the right-hand side
    ax.text(1.03, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

#################
# MAIN FUNCTION #
#################

def main():
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(duration_result_dir, f'*/*RUN{num_run}*.npz')))

    # get used DUs
    du_list = get_npz_dus(file_list)

    # loop through used DUs
    for du in du_list:
        get_duration(du)

if __name__ == "__main__":
    with record_run_time():
        main()