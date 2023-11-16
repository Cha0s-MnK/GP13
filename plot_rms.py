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

# loop through all good DUs to load the NPZ file and plot the histograms
for du in du_list:
    # load corresponding NPZ file
    rms_name     = f'DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_cutoff{cutoff_frequency}.npz'
    rms_npz_file = os.path.join(rms_result_dir, rms_name)
    data         = np.load(rms_npz_file, allow_pickle=True)
    
    # create a figure and adjust figure size
    plt.figure(figsize=(36, 12), dpi=256)
    plt.suptitle(f'Noise Distribution of DU{du}')
    
    # plot a histogram for each channel in subfigures
    for id, channel in enumerate(channels, 1):
        plt.subplot(3, 1, id)  # 3 rows, 1 column layout
        plt.hist(data[f'noise_list{channel}'], bins=50, alpha=0.75, label=f'Channel {channel}')
        plt.xlabel('Noise value')
        plt.ylabel('Frequency')
        plt.title(f'Channel {channel}', fontsize=20)
        plt.legend()
    
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