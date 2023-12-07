#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

name = f'{trace_wanted}-rms-psd'

# create the save directory if it does not exist
save_dir = os.path.join(plot_dir, name)
os.makedirs(save_dir, exist_ok=True)

# define a dictionary mapping DUs to colours and line styles
du_styles = {
    1010: {'colour': 'green', 'linestyle': '-', 'offset': 110},
    1013: {'colour': 'yellow', 'linestyle': '-', 'offset': 100},
    1017: {'colour': 'pink', 'linestyle': '-', 'offset': 90},
    1019: {'colour': 'black', 'linestyle': '-', 'offset': 80},
    1020: {'colour': 'orange', 'linestyle': '-', 'offset': 70},
    1021: {'colour': 'purple', 'linestyle': '-', 'offset': 60},
    1029: {'colour': 'cyan', 'linestyle': '-', 'offset': 50},
    1032: {'colour': 'brown', 'linestyle': '-', 'offset': 40},
    1035: {'colour': 'magenta', 'linestyle': '-', 'offset': 30},
    1041: {'colour': 'gray', 'linestyle': '-', 'offset': 20},
    1076: {'colour': 'red', 'linestyle': '-', 'offset': 0},
    1085: {'colour': 'blue', 'linestyle': '-', 'offset': 0}
}

# adjust the figure size
fig_length = 24
fig_width  = 12
fig_dpi    = 165

##################
# CORE FUNCTIONS #
##################

# plot RMS / mean PSDs of different DUs and channels for 1 day
def supplot1day_rms_psd(date, du_list):
    # update NPZ files according to current date
    file_list = get_npz_files(name=f'{trace_wanted}-rms-psd', date=date)

    # make dictionaries for easier indexing
    rmses_list     = {du: {} for du in du_list}
    mean_psds_list = {du: {} for du in du_list}
    times_list     = {du: {} for du in du_list}

    # loop through used DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            rmses_list[du][channel]     = []
            mean_psds_list[du][channel] = []
            times_list[du][channel]     = []

    # loop through NPZ files
    for file in file_list:
        # get information from this NPZ file
        npz_file = np.load(file, allow_pickle=True)
        du       = get_npz_du(file)
        if du not in du_list:
            continue
        for channel in channels:
            rmses_list[du][channel]     = npz_file[f'rmses_list{channel}']
            mean_psds_list[du][channel] = npz_file[f'mean_psds_list{channel}']
            times_list[du][channel]     = npz_file[f'times_list{channel}']

    plot1day_rms_psd(date=date, du_list=du_list, datasetX=times_list, datasetY=rmses_list, labelY='RMS of the trace / ADC units', filename='rms-hour')
    plot1day_rms_psd(date=date, du_list=du_list, datasetX=times_list, datasetY=mean_psds_list, labelY='Mean PSDs of the trace from 112.3MHz to 146.2MHz / $V^2 MHz^{-1}$', filename='mean-psd-hour')

def plot1day_rms_psd(date, du_list, datasetX, datasetY, labelY, filename):
    # create a layout for the subplots
    fig, axes = plt.subplots(3, 1, figsize=(fig_length, fig_width), dpi=fig_dpi)
    
    # Plot the transient/pulse rates of the filtered trace per channel
    day_start = datetime.combine(datetime.strptime(date, '%Y%m%d').date(), time())
    for channel_id, channel in enumerate(channels):
        subplot1channel_rms_psd(ax=axes[channel_id], channel=channel, datasetX=datasetX, datasetY=datasetY, day_start=day_start, du_list=du_list)

    # set common labels and title
    fig.text(0.5, 0.01, 'CST', ha='center', fontsize=18)
    fig.text(0.01, 0.5, labelY, va='center', rotation='vertical', fontsize=18)
    fig.suptitle(f'{date}', fontsize=20)

    # adjust the layout
    plt.tight_layout(rect=[0.02, 0.02, 0.99, 0.99])

    # Save the figure as a PNG file
    dus = '-'.join(str(du) for du in du_list)
    save_file = os.path.join(save_dir, f'{trace_wanted}-{filename}_RUN{nums_run}_DU{dus}_{date}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')

def subplot1channel_rms_psd(ax, channel, datasetX, datasetY, day_start, du_list):
    # plot used DUs for this channel
    for du in du_list:
        # get styles of current DU
        style       = du_styles.get(du)
        colour      = style['colour']
        linestyle   = style['linestyle']
        offset      = style['offset']

        # convert lists to arrays
        dataX = np.array(datasetX[du][channel])
        dataY = np.array(datasetY[du][channel]) + offset

        ax.plot(dataX, dataY, marker='o', color=colour, linestyle=linestyle, label=f'DU{du} channel {channel}, offset = {offset}')

    # set X-ticks
    #ax.set_xlim(day_start, day_start + timedelta(days=1))
    ax.set_xlim(day_start + timedelta(hours=8), day_start + timedelta(hours=9))
    #ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 10, 20, 30, 40, 50]))
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    # set Y-scale
    ax.set_yscale('log')

    # set the title on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(ncol=3, frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

def plot1day_psd_band(date, du_list):
    # load data to plot
    rmses_list, psd_means_list, psd_bands_list, times_list = load1day_rms_psd(date=date, du_list=du_list)

    # choose the start and stop time
    day_start  = datetime.combine(datetime.strptime(date, '%Y%m%d'), time())
    start_time = day_start + timedelta(hours=8)
    stop_time  = day_start + timedelta(hours=9)

    # plot a PNG file for each DU
    for du in du_list:
        fig, axes = plt.subplots(3, 1, figsize=(fig_length, fig_width), dpi=fig_dpi)
    
        for channel_id, channel in enumerate(channels):
            # use a mask to filter the data
            time_mask            = (times_list[du][channel] >= start_time) & (times_list[du][channel] < stop_time)
            filter_time_list     = times_list[du][channel][time_mask]
            filter_rms_list      = np.array(rmses_list[du][channel])[time_mask]
            filter_psd_mean_list = np.array(psd_means_list[du][channel])[time_mask]
            filter_psd_band_list = np.array(psd_bands_list[du][channel])[time_mask]

            subplot1day_psd_band(ax=axes[channel_id], channel=channel, dataX=filter_time_list, dataZ=filter_psd_band_list)

        # add a color bar
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        sm = plt.cm.ScalarMappable(cmap='RdYlBu', norm=LogNorm(vmin=1e-14, vmax=1e-6))
        sm.set_array([])
        fig.colorbar(sm, cax=cbar_ax)

        # set common labels and title
        fig.text(0.5, 0.01, 'CST', ha='center', fontsize=18)
        fig.text(0.01, 0.5, 'PSDs of the trace (115-144MHz) / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
        fig.suptitle(f'DU{du} on {date}', fontsize=20)

        # adjust the layout: left, bottom, right, top
        plt.tight_layout(rect=[0.02, 0.02, 0.99, 0.99])

        # save as a PNG file
        save_file = os.path.join(save_dir, f'{wanted_trace}-psd-band_RUN{nums_run}_DU{du}_{date}.png')
        plt.savefig(save_file)
        print(f'Saved: {save_file}')

def subplot1day_psd_band(ax, channel, dataX, dataZ):
    dataY = fft_frequency[470:590]

    # plot the 2-D colour map
    gridX, gridY = np.meshgrid(dataX, dataY)
    gridZ        = dataZ.T
    #colour_map2d = ax.pcolormesh(gridX, gridY, gridZ, shading='auto', norm=LogNorm(), cmap='jet')
    colour_map2d = ax.pcolormesh(gridX, gridY, gridZ, shading='auto', norm=LogNorm(vmin=1e-14, vmax=1e-6), cmap='RdYlBu')

    # set X-ticks
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    # set the channel on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

def load1day_rms_psd(date, du_list):
    file_list = get_npz_files(name=name, date_list=[date], du_list=du_list)

    # make dictionaries for easier indexing and initiate them
    rmses_list     = {du: {} for du in du_list}
    psd_means_list = {du: {} for du in du_list}
    psd_bands_list = {du: {} for du in du_list}
    times_list     = {du: {} for du in du_list}
    for du in du_list:
        for channel in channels:
            rmses_list[du][channel]     = []
            psd_means_list[du][channel] = []
            psd_bands_list[du][channel] = []
            times_list[du][channel]     = []

    # get info from NPZ files
    for file in file_list:
        npz_file = np.load(file, allow_pickle=True)
        du       = get_npz_du(file)
        for channel in channels:
            rmses_list[du][channel]     = npz_file[f'rmses_list{channel}']
            psd_means_list[du][channel] = npz_file[f'psd_means_list{channel}']
            psd_bands_list[du][channel] = npz_file[f'psd_bands_list{channel}']
            times_list[du][channel]     = npz_file[f'times_list{channel}']
    
    return rmses_list, psd_means_list, psd_bands_list, times_list

#################
# MAIN FUNCTION #
#################

def main():
    #file_list = get_npz_files(name=name)
    date_list = get_npz_dates(file_list=get_npz_files(name=name, date_list=False))
    du_list   = get_npz_dus(file_list=get_npz_files(name=name, du_list=False))

    for date in date_list:
        plot1day_psd_band(date=date, du_list=du_list)

if __name__ == "__main__":
    with record_run_time():
        main()