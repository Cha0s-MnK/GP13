#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

# create the save directory if it does not exist
save_dir = os.path.join(plot_dir, 'mean-psd')
os.makedirs(save_dir, exist_ok=True)

# define a dictionary mapping DUs to colours and line styles
du_styles = {
    1010: {'colour': 'green', 'linestyle': '-', 'offset': 120},
    1013: {'colour': 'yellow', 'linestyle': '--', 'offset': 110},
    1017: {'colour': 'pink', 'linestyle': '-.', 'offset': 100},
    1019: {'colour': 'black', 'linestyle': ':', 'offset': 90},
    1020: {'colour': 'orange', 'linestyle': '-', 'offset': 80},
    1021: {'colour': 'purple', 'linestyle': '--', 'offset': 70},
    1029: {'colour': 'cyan', 'linestyle': '-.', 'offset': 60},
    1031: {'colour': 'darkorange', 'linestyle': ':', 'offset': 50},
    1032: {'colour': 'brown', 'linestyle': '-', 'offset': 40},
    1035: {'colour': 'magenta', 'linestyle': '--', 'offset': 30},
    1041: {'colour': 'gray', 'linestyle': '-.', 'offset': 20},
    1075: {'colour': 'red', 'linestyle': ':', 'offset': 10},
    1085: {'colour': 'blue', 'linestyle': '-', 'offset': 0}
}

##################
# CORE FUNCTIONS #
##################

def plot1day_mean_psd(date, du_list):
    # load data
    start_time = datetime.combine(datetime.strptime(date, '%Y%m%d'), time()) + timedelta(hours=8)
    stop_time  = start_time + timedelta(days=1)
    file_list  = get_npz_files(name='rms-psd', date_list=[date], du_list=du_list)
    du_list    = get_npz_dus(file_list=file_list)
    times_list, rmses_list, psds_list, temperatures_list, mean_psds = load1day_rms_psd(file_list=file_list, du_list=du_list, start_time=start_time, stop_time=stop_time)

    fig, axes = plt.subplots(3, 1, figsize=(16, 24), dpi=165)

    for channel_id, channel in enumerate(channels):
        subplot1channel_mean_psd(ax=axes[channel_id], channel=channel, du_list=du_list, datasetY=mean_psds)

    # set common labels and title
    fig.text(0.5, 0.01, 'Frequency / MHz', ha='center', fontsize=18)
    fig.text(0.01, 0.5, 'PSDs (50-200MHz) / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
    #fig.suptitle(f'DU{du} on {date}', fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.02, 0.02, 0.95, 1.0])

    # save as a PNG file
    save_file = os.path.join(save_dir, f'mean-psd_RUN{nums_run}_{len(du_list)}DUs_{date}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')

def subplot1channel_mean_psd(ax, channel, du_list, datasetY):
    for du in du_list:
        # get styles of current DU
        style      = du_styles.get(du)
        colour     = style['colour']
        linestyle  = style['linestyle']
        offset     = style['offset']

        ax.plot(fft_frequency[205:820], datasetY[du][channel], color=colour, linestyle=linestyle, label=f'DU{du}')

    # load galactic noise simulation for this channel
    galaxy_sim = np.load(os.path.join(galaxy_dir, galaxy_name[channel_mapping[channel]]))

    # apply noise according to the hour and correct the 20dB linear gain
    galaxy_noise = galaxy_sim[:,2] / linear_gain / linear_gain

    # plot the galactic noise simulation
    ax.semilogy(range(10, 250), galaxy_noise[0:240], color='black', linestyle='--', label='Galactic Noise')

    # add vertical lines at 30MHz and 50MHz
    #ax.axvline(x=30, color='purple', linestyle=':', linewidth=2, label='30 MHz')
    #ax.axvline(x=50, color='green', linestyle=':', linewidth=2, label='50 MHz')

    # set X-ticks
    ax.set_xlim([50, 200])
    ax.set_xticks(np.arange(50, 201, 10))

    ax.set_yscale('log')
    ax.set_ylim([1e-12, 1e-7])

    # set the title on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    ax.grid(True)
    ax.legend(ncol=3, frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

#################
# MAIN FUNCTION #
#################

def main():
    date_list = get_npz_dates(file_list=get_npz_files(name='rms-psd', date_list=False))

    for date in date_list:
        plot1day_mean_psd(date=date, du_list=[1010, 1013, 1017, 1019, 1020, 1021, 1029, 1032, 1035, 1041, 1085])

if __name__ == "__main__":
    with record_run_time():
        main()