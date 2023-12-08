#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

# create the save directory if it does not exist
save_dir = os.path.join(plot_dir, f'{trace_wanted}-rms-psd')
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

def plot1day_linear(date, du_list):
    # update NPZ files according to current date
    file_list = get_npz_files(name=f'{trace_wanted}-rms-psd', date=date)

    # make dictionaries for easier indexing
    rmses_list       = {du: {} for du in du_list}
    mean_psds_list   = {du: {} for du in du_list}
    times_list       = {du: {} for du in du_list}
    common_times_set = {channel: None for channel in channels}

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

    for channel in channels:
        for du in du_list:
            if common_times_set[channel] is None:
                common_times_set[channel] = set(times_list[du][channel])
            else:
                common_times_set[channel].intersection_update(times_list[du][channel])
        
        # convert the set to a sorted list
        common_times_list = sorted(list(common_times_set[channel]))
        
        for du in du_list:
            # find the indices of the common times in the original time list
            indices = np.isin(times_list[du][channel], common_times_list)
        
            # filter using the indices
            rmses_list[du][channel] = rmses_list[du][channel][indices]
            mean_psds_list[du][channel] = mean_psds_list[du][channel][indices]
    
    # create a layout for the subplots and flatten the axes array for easier indexing
    fig, ax = plt.subplots(figsize=(10, 10), dpi=165)

    # plot the transient/pulse rates of the filtered trace per channel
    subplot1channel_linear(ax=ax, dataX=rmses_list[1076]['X'], dataY=mean_psds_list[1085]['X'])

    # set common labels and title
    fig.text(0.5, 0.01, 'RMS of traces from DU1076 / ADC units', ha='center', fontsize=18)
    #fig.text(0.5, 0.01, 'Mean PSD of traces from DU1076 / $V^2 MHz^{-1}$', ha='center', fontsize=18)
    #fig.text(0.01, 0.5, 'RMS of traces from DU1085 / ADC units', va='center', rotation='vertical', fontsize=18)
    fig.text(0.01, 0.5, 'Mean PSD of traces from DU1085 / $V^2 MHz^{-1}$', va='center', rotation='vertical', fontsize=18)
    fig.suptitle(f'Comparison in channel X on {date}', fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.02, 0.02, 0.99, 0.99])

    # save the figure as a PNG file
    save_file = os.path.join(save_dir, f'{trace_wanted}-linear_RUN{nums_run}_{date}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')

def subplot1channel_linear(ax, dataX, dataY):
    # take the log of the original data
    dataX = np.log(np.array(dataX))
    dataY = np.log(np.array(dataY))

    # plot used DUs for this channel
    ax.scatter(dataX, dataY, color='blue', marker='o', label=f'Original Data')

    # perform linear regression and create a function using the slope and intercept
    slope, intercept = np.polyfit(dataX, dataY, 1)
    fit_line = slope * dataX + intercept
    ax.plot(dataX, fit_line, color='red', linestyle='-', label=f'Fit Line (y = {slope:.2f}x + {intercept:.2f})')

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

#################
# MAIN FUNCTION #
#################

def main():
    # get NPZ files
    file_list = get_npz_files(name=f'{trace_wanted}-rms-psd')

    # get the dates
    date_list = get_npzs_dates(file_list=file_list)
    date_list = ['20231122', '20231123', '20231124']

    # input plotted DUs manually
    du_list = [1076, 1085]

    # loop through the dates
    for date in date_list:
        # plot for 1 day
        plot1day_linear(date=date, du_list=du_list)

if __name__ == "__main__":
    with record_run_time():
        main()