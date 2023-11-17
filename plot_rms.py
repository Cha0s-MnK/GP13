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
    du = int(basename.split('_')[0][2:])
    du_set.add(du)

# convert the set to a list in order
du_list = sorted(list(du_set))

# only consider good DUs
du_list = [du for du in du_list if du in good_du_list]

# print the list of all used DUs
print(f'\nNPZ files contain data from following DUs: \n{du_list}\n')

#######################
# PLOT RMS TO ANALYZE #
#######################

# loop through all good DUs to load the NPZ file and plot the histograms
for du in du_list:
    # load corresponding NPZ file
    rms_name     = f'DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_cutoff{cutoff_frequency}.npz'
    rms_npz_file = os.path.join(rms_result_dir, rms_name)
    data         = np.load(rms_npz_file, allow_pickle=True)
    
    # create a figure and adjust figure size
    fig = plt.figure(figsize=(12, 20), dpi=256)
    fig.suptitle(f'Noise Distribution of DU{du}', fontsize=20)
    
    # plot a histogram for each channel in subfigures
    for id, channel in enumerate(channels, 1):
        # create a 3 rows, 1 column layout
        ax = plt.subplot(3, 1, id)

        # extract the noise data for the channel
        noise_channel = data[f'noise_list{channel}']

        # plot the histogram and calculate the parameters
        count, bins, ignored = ax.hist(noise_channel, bins=100, density=True, color='blue', alpha=0.75, label=f'Channel {channel}')
    
        # plot a Gaussian distribution compared with the noise data
        sigma = np.std(noise_channel)
        normal_distribution = norm.pdf(bins, 0, sigma)
        ax.plot(bins, normal_distribution, linewidth=2, color='r', label=f'Gaussian (sigma={sigma:.2f})')

        # determine the range for X-ticks and set the X-ticks symmetrically around 0
        max_abs_value = max(abs(noise_channel.max()), abs(noise_channel.min()))
        Xtick_step = max_abs_value / 6
        Xticks = np.arange(-max_abs_value, max_abs_value + Xtick_step, Xtick_step)
        ax.set_xticks(Xticks)

        # enable the grid
        ax.grid(True)

        # set individual title on the right-hand side
        ax.text(1.01, 0.5, f'Channel {channel}', verticalalignment='center', horizontalalignment='left', transform=ax.transAxes, fontsize=18, rotation=-90)

        # enable the legend
        ax.legend()
    
    # add common X and Y labels
    fig.text(0.5, 0.04, 'Noise Value / ADC Counts', ha='center', va='center', fontsize=16)
    fig.text(0.04, 0.5, 'Probability Density', ha='center', va='center', rotation='vertical', fontsize=16)

    # adjust layout: left, bottom, right, top
    plt.tight_layout(rect=[0.04, 0.04, 1.0, 1.0]) 

    # save the figure as a PNG file
    rms_png_name = f'rms_DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_cutoff{cutoff_frequency}.png'
    rms_png_file = os.path.join(rms_plot_dir, rms_png_name)
    plt.savefig(rms_png_file)
    print(f'Saved: {rms_png_file}')

    # close the figure to free up memory
    plt.close()

# record the running time
end_time = wall_time.perf_counter()
run_time = end_time - start_time
print(f"\nWhole program executed in: {run_time} seconds")