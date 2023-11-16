#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#############################
# GET NPZ FILES AND DU LIST #
#############################

# get selected NPZ files
file_list = sorted(glob.glob(rms_npz_files))

# create an empty set to store unique DUs
du_set = set()

# show all NPZ files and extract DUs
print('\nData from following NPZ files:')
for file in file_list:
    basename = os.path.basename(file)
    print(basename)
    
    # split the filename on underscore
    du = basename.split('_')[0][2:]
    du_set.add(du)

# convert the set to a list in order
du_list = sorted(list(du_set))

# print the list of all used DUs
print(f'\nNPZ files contain data from following DUs: \n{du_list}\n')

#######################
# PLOT RMS TO ANALYZE #
#######################

# loop through all DUs to plot transient rates for 3 ADC channels
for du in du_list:
    # create a 3x1 layout for the subplots and adjust figure size
    fig, axes = plt.subplots(3, 1, figsize=(36, 12), dpi=256)

    # flatten axes array for easier indexing
    axes = axes.flatten()

    # plot the data per channel
    for id, channel in enumerate(channels):
        # reformat DunHuang hours and transient rates to plot
        hours_list[channel][du] = list(itertools.chain(*hours_list[channel][du]))
        rates_list[channel][du] = 1+np.array(list(itertools.chain(*rates_list[channel][du])))

        # locate corresponding axis
        ax = axes[id]

        # scatter points and add vertical lines to make the plot more clearly
        ax.scatter(hours_list[channel][du], rates_list[channel][du], color='blue', s=9, alpha=0.6)
        ax.vlines(hours_list[channel][du], 0, rates_list[channel][du], color='blue', linestyle='solid', linewidth=4, alpha=0.75)

        # set axis scales
        ax.set_yscale('log')

        # set only the lower limit of the y-axis
        ax.set_ylim(bottom=1)

        # enable gird
        ax.grid(True)

        # enable ticks on both left and right sides of the plot
        ax.tick_params(axis='y', labelleft=True, labelright=True)

        # set title on the right-hand side
        ax.text(1.01, 0.5, f'Channel {channel}', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=20, rotation=-90)
        
        # set X-limits
        ax.set_xlim(datetime(2023, 10, 11, 6), datetime(2023, 10, 31, 18))

        # set the major tick locator to every 12 hours
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))

        # only set X-tick labels for the last subplot
        if id < len(channels) - 1:
            ax.set_xticklabels([])
        else:
            # set the formatting of the ticks to include the hour
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))

            # rotate X labels
            for tick in ax.get_xticklabels():
                tick.set_rotation(15)

    # set X label and Y label for entire figure
    fig.text(0.5, 0.0, 'DunHuang Time in Date and Hour', ha='center', fontsize=20)
    fig.text(0.0, 0.5, 'Transient Rate / kHz', va='center', rotation='vertical', fontsize=20)
    
    # add the figure legend
    fig.legend(custom_lines, ['Sunrise', 'Sunset'], loc='upper right', fontsize=18, bbox_to_anchor=(0.98,1), bbox_transform=plt.gcf().transFigure)

    # add a main/super title for the entire figure
    plt.suptitle(f'Time Evolution of Transient Rates \n DU{du}, Threshold = {num_threshold}, Separation = {standard_separation}, Crossing = {num_crossings}, Cut-off frequency = {cutoff_frequency}, Noises = [{noises[0]}, {noises[1]}, {noises[2]}]', fontsize=24)

    # adjust layout
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 1.0])

    # save the figure as a PNG file
    rms_png_name = f'rate_DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_cutoff{cutoff_frequency}.png'
    rms_png_file = os.path.join(rms_png_dir, rms_png_name)
    plt.savefig(os.path.join(rate_dir, rate_file))
    print(f'Saved: {rate_dir}{rate_file}')

    # close the figure to free up memory
    plt.close(fig)

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")