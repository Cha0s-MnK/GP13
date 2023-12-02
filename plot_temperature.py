#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

# use global variables
du     = 0
subdir = 'temperature'

##################
# PLOT FUNCTIONS #
##################

def plot1du_temperature():
    # update NPZ files according to current DU
    file_list = sorted(glob(os.path.join(result_dir, subdir, f'normal-temperature*RUN{num_run}_DU{du}*.npz')))

    # make dictionaries for easier indexing and initiate them with empty lists
    noises_list              = {channel: [] for channel in channels}
    temperature_list         = []
    filter_noises_list       = {channel: [] for channel in channels}
    filter_temperatures_list = {channel: [] for channel in channels}

    # loop through NPZ files
    for file in file_list:
        # get information from this file
        npz_file = np.load(file, allow_pickle=True)
        for channel in channels:
            test = npz_file[f'noises_list{channel}']
            noises_list[channel].extend(npz_file[f'noises_list{channel}'])
        temperature_list.extend(npz_file['temperature_list'])
    
    # filter abnormal background noises
    for channel in channels:
        filter_temperatures_list[channel], filter_noises_list[channel] = abnormal_noise_filter(temperature_list=temperature_list, noise_list=noises_list[channel])

    my_plot(row=3, column=1, datasetX=filter_temperatures_list, datasetY=filter_noises_list, labelX='Temperature', labelY='Log of Background Noise / ADC units', title=f'Noise-Temperature Curve for DU{du}')

def my_plot(row, column, datasetX, datasetY, labelX, labelY, title):
    # create a layout for the subplots
    fig, axes = plt.subplots(row, column, figsize=(12, 24), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the transient/pulse rates of the filtered trace per channel
    for id, channel in enumerate(channels):
        my_subplot(axes[id], channel, datasetX[channel], datasetY[channel])

    # set common labels and title
    fig.text(0.5, 0.0, labelX, ha='center', fontsize=18)
    fig.text(0.0, 0.5, labelY, va='center', rotation='vertical', fontsize=18)
    plt.suptitle(title, fontsize=20)

    # adjust the layout: left, bottom, right, top
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])

    # create the save directory if it does not exist
    save_dir = os.path.join(plot_dir, subdir)
    os.makedirs(save_dir, exist_ok=True)

    # save the figure as a PNG file
    save_file = os.path.join(save_dir, f'normal-temperature_RUN{num_run}_DU{du}.png')
    plt.savefig(save_file)
    print(f'Saved: {save_file}')

    # close the figure to free up memory
    plt.close(fig)

def my_subplot(ax, channel, dataX, dataY):
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
    #ax.set_yscale('log')
    #ax.set_ylim([5e-2, 1e6])

    # set the title on the right-hand side
    ax.text(1.03, 0.5, f'channel {channel} ', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=16, rotation=-90)

    # enable grid, legend and Y-ticks on both sides
    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', fontsize=14)
    ax.tick_params(axis='y', labelleft=True, labelright=True)

def abnormal_noise_filter(temperature_list, noise_list):
    '''
    # convert lists to arrays
    temperature_array = np.array(temperature_list)
    noise_array       = np.array(noise_list)

    mean = np.mean(noise_array)
    mask = noise_array < 2 * mean

    # mask the arrays
    temperature_list = list(noise_array[mask])
    noise_list       = list(temperature_array[mask])
    '''
    return temperature_list, noise_list

######################
# MAIN PLOT FUNCTION #
######################

def main():
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(result_dir, subdir, f'normal-temperature*RUN{num_run}*.npz')))

    # get used DUs
    du_list = get_npz_dus(file_list)

    # declare global variables
    global du

    # loop through used DUs
    for cur_du in du_list:
        du = cur_du
        plot1du_temperature()

if __name__ == "__main__":
    with record_run_time():
        main()