#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

# create the save directory if it does not exist
name = 'rms-psd'
save_dir = os.path.join(plot_dir, name)
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

# adjust the figure size
fig_length = 150
fig_width  = 24
fig_dpi    = 165

##################
# CORE FUNCTIONS #
##################

def plot1day_rms(date, du_list):
    # load data
    start_time = datetime.combine(datetime.strptime(date, '%Y%m%d'), time()) + timedelta(hours=8) # convert UTC to CST
    stop_time  = start_time + timedelta(days=1)
    file_list  = get_npz_files(name=name, date_list=[date], du_list=du_list)
    du_list    = get_npz_dus(file_list=file_list)
    times_list, rmses_list, psds_list, temperatures_list, psd_means_list = load1day_rms_psd(file_list=file_list, du_list=du_list, start_time=start_time, stop_time=stop_time)
    
    fig, axes = plt.subplots(3, 1, figsize=(fig_length, fig_width), dpi=fig_dpi)
    
    for channel_id, channel in enumerate(channels):
        subplot1channel_rms(ax=axes[channel_id], channel=channel, du_list=du_list, datasetX=times_list, datasetY=rmses_list, start_time=start_time, stop_time=stop_time)

    # set common labels and title
    fig.text(0.5, 0.01, 'CST', ha='center', fontsize=18)
    fig.text(0.01, 0.5, 'RMS of the trace / ADC units', va='center', rotation='vertical', fontsize=18)
    fig.suptitle(f'{date}', fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.02, 0.02, 0.99, 0.99])

    # save as a PNG file
    save_file = os.path.join(save_dir, f'rms_RUN{nums_run}_{len(du_list)}DUs_{date}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')

def subplot1channel_rms(ax, channel, du_list, datasetX, datasetY, start_time, stop_time):
    for du in du_list:
        # get styles of current DU
        style      = du_styles.get(du)
        colour     = style['colour']
        linestyle  = style['linestyle']
        offset     = style['offset']

        # convert lists to arrays
        dataX = np.array(datasetX[du][channel])
        dataY = np.array(datasetY[du][channel]) + offset

        ax.plot(dataX, dataY, marker='o', color=colour, linestyle=linestyle, label=f'DU{du} channel {channel}, offset = {offset}')

    # set X-ticks
    ax.set_xlim([start_time, stop_time])
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 10, 20, 30, 40, 50]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ax.set_ylim([0, 180])

    # set the title on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(ncol=3, frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

def plot1day_psd(date, du_list):
    # load data
    start_time = datetime.combine(datetime.strptime(date, '%Y%m%d'), time()) + timedelta(hours=8) + timedelta(hours=8)
    stop_time  = start_time + timedelta(hours=3)
    file_list  = get_npz_files(name=name, date_list=[date], du_list=du_list)
    du_list    = get_npz_dus(file_list=file_list)
    times_list, rmses_list, psds_list, temperatures_list, psd_means_list = load1day_rms_psd(file_list=file_list, du_list=du_list, start_time=start_time, stop_time=stop_time)

    # plot a PNG file for each DU
    for du in du_list:
        fig, axes = plt.subplots(3, 1, figsize=(64, 48), dpi=165)
    
        for channel_id, channel in enumerate(channels):
            subplot1day_psd_band(ax=axes[channel_id], channel=channel, dataX=times_list[du][channel], dataZ=psds_list[du][channel])

        # add a color bar
        cbar_ax = fig.add_axes([0.95, 0.05, 0.02, 0.9]) # left, bottom, width, height
        sm = plt.cm.ScalarMappable(cmap='jet', norm=LogNorm(vmin=1e-14, vmax=1e-6))
        sm.set_array([])
        fig.colorbar(sm, cax=cbar_ax)

        # set common labels and title
        fig.text(0.5, 0.01, 'CST', ha='center', fontsize=18)
        fig.text(0.01, 0.5, 'PSDs of the trace (50-200MHz) / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
        #fig.suptitle(f'DU{du} on {date}', fontsize=20)

        # adjust the layout: left, bottom, right, top
        plt.tight_layout(rect=[0.02, 0.02, 0.95, 1.0])

        # save as a PNG file
        save_file = os.path.join(save_dir, f'psd_RUN{nums_run}_DU{du}_{date}.png')
        plt.savefig(save_file)
        print(f'Saved: {save_file}')

def subplot1day_psd_band(ax, channel, dataX, dataZ):
    # plot the 2-D colour map
    dataY        = fft_frequency[205:820]
    gridX, gridY = np.meshgrid(np.array(dataX), np.array(dataY))
    gridZ        = np.array(dataZ).T
    colour_map2d = ax.pcolormesh(gridX, gridY, gridZ, shading='auto', norm=LogNorm(vmin=1e-14, vmax=1e-6), cmap='jet')

    # set X-ticks
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    # set the channel on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    ax.grid(True)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

def plot1day_temperature(date, du_list):
    # load data
    times_list, rmses_list, psd_means1_list, psd_means2_list, psd_means3_list, psd_bands_list, temperatures_list = load1day_rms_psd(date=date, du_list=du_list)

    # plot a PNG file for each DU
    for du in du_list:
        fig, axes = plt.subplots(3, 1, figsize=(fig_length, fig_width), dpi=fig_dpi)

        for channel_id, channel in enumerate(channels):
            subplot1day_temperature(ax=axes[channel_id], channel=channel, dataX=temperature_array, dataY=psd_mean1_array)

        # set common labels and title
        fig.text(0.5, 0.01, 'Temperature / C', ha='center', fontsize=18)
        fig.text(0.01, 0.5, 'Log of Mean PSD / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
        fig.suptitle(f'DU{du} on {date}', fontsize=20)

        # adjust the layout: left, bottom, right, top
        plt.tight_layout(rect=[0.02, 0.02, 0.95, 0.99])

        # save as a PNG file
        save_file = os.path.join(save_dir, f'{wanted}-temperature_RUN{nums_run}_DU{du}_{date}.png')
        plt.savefig(save_file)
        print(f'Saved: {save_file}')

def subplot1channel_temperature(ax, channel, datasetX, datasetY, day_start, du_list):
    # convert lists to arrays
    dataX = np.array(dataX)
    dataY = np.log(np.array(dataY))

    # perform linear regression and create a function using the slope and intercept
    slope, intercept = np.polyfit(dataX, dataY, 1)
    fit_line = slope * dataX + intercept

    # plot the data and the fit line
    ax.scatter(dataX, dataY, color='blue', label='Original Data')
    ax.plot(dataX, fit_line, color='red', label=f'Fit Line (y = {slope:.2f}x + {intercept:.2f})')

    # set Y-scale
    ax.set_yscale('log')

    # set the title on the right-hand side
    ax.text(1.02, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

def load1day_rms_psd(file_list, du_list, start_time, stop_time):
    # make dictionaries for easier indexing and initiate them
    rmses_list        = {du: {} for du in du_list}
    psds_list         = {du: {} for du in du_list}
    times_list        = {du: {} for du in du_list}
    temperatures_list = {du: {} for du in du_list}
    psd_means_list    = {du: {} for du in du_list}
    for du in du_list:
        for channel in channels:
            rmses_list[du][channel]        = []
            psds_list[du][channel]         = []
            times_list[du][channel]        = []
            temperatures_list[du][channel] = []
            psd_means_list[du][channel]    = []

    # get info from NPZ files
    for file in file_list:
        npz_file = np.load(file, allow_pickle=True)
        du       = get_npz_du(file)
        for channel in channels:
            rmses_list[du][channel]        = npz_file[f'rmses_list{channel}']
            psds_list[du][channel]         = npz_file[f'psds_list{channel}']
            times_list[du][channel]        = npz_file[f'times_list{channel}']
            temperatures_list[du][channel] = npz_file[f'temperatures_list{channel}']

            psd_means_list[du][channel]    = [np.mean(sublist) for sublist in psds_list[du][channel]]

    # use a mask to filter the data
    for file in file_list:
        for du in du_list:
            for channel in channels:
                time_mask                      = (times_list[du][channel] >= start_time) & (times_list[du][channel] < stop_time)
                times_list[du][channel]        = times_list[du][channel][time_mask]
                rmses_list[du][channel]        = rmses_list[du][channel][time_mask]
                psds_list[du][channel]         = psds_list[du][channel][time_mask]
                temperatures_list[du][channel] = temperatures_list[du][channel][time_mask]
                psd_means_list[du][channel]    = np.array(psd_means_list[du][channel])[time_mask]
    
    return times_list, rmses_list, psds_list, temperatures_list, psd_means_list

#################
# MAIN FUNCTION #
#################

def main():
    #file_list = get_npz_files(name=name)
    date_list = get_npz_dates(file_list=get_npz_files(name=name, date_list=False))

    for date in date_list:
        plot1day_psd(date=date, du_list=False)

if __name__ == "__main__":
    with record_run_time():
        main()