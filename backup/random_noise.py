#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

##################################
# SIMULATION OF A GAUSSIAN NOISE #
##################################

# generate a Gaussian noise
np.random.seed(0) # for reproducibility
noise_trace = np.random.normal(loc=0.0, scale=1.0, size=num_samples)

# calculate the threshold and search for time windows
threshold = num_threshold * np.std(noise_trace)
window_list = search_windows(noise_trace, threshold, standard_separation)

# create a figure with a specific size and plot the noise trace
fig, ax = plt.subplots(figsize=(36, 4), dpi=256)
ax.plot(time_axis, noise_trace, color='blue', label='noise trace')

# add horizontal dashed lines at +/-Threshold
ax.axhline(y=threshold, color='red', linestyle='--', label=f'+/-Threshold={threshold:.2f}')
ax.axhline(y=-threshold, color='red', linestyle='--')

# set X-limits
ax.set_xlim(-4, 4100)

# highlight the time windows in the time axis
for i, window in enumerate(window_list):
    start_id, stop_id = window # 1 sample = 2 nanoseconds
    start_time = 2 * start_id # 1 sample = 2 nanoseconds
    stop_time  = 2 * stop_id
    # label only the 1st axvspan
    if i == 0:
        ax.axvspan(start_time, stop_time, color='green', alpha=0.3, label='Time Window(s)')
    else:
        ax.axvspan(start_time, stop_time, color='green', alpha=0.3)

    # annotate time windows
    centre_time = (start_time + stop_time) / 2
    ylim = ax.get_ylim()
    ax.text(centre_time, 1.02*ylim[1], f'[{start_time:.1f}, {stop_time:.1f} ]', color='black', fontsize=10, ha='center', va='bottom')

# adjust the size of the x-tick and y-tick labels
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)

# set an x-tick every 250 units
ax.set_xticks(np.arange(min(time_axis), max(time_axis)+1, 250))

# add some comments
ax.legend(loc='upper right')
ax.set_xlabel('Time / ns', fontsize=16)
ax.set_ylabel('ADC Counts', fontsize=16)
ax.set_title(f'Gaussian Noise Trace with Detected Windows \n Threshold = {num_threshold}, Separation = {standard_separation}', fontsize=20, pad=20)

# adjust layout
plt.tight_layout()

# save the figure as a PNG file
plt.savefig('plot/random_noise.png')

######################################
# SIMULATION OF DIFFERENT THRESHOLDS #
######################################

# 1 simulation of a noise trace
def noise_sim(num_threshold):
    # generate a Gaussian noise trace
    noise_trace = np.random.normal(loc=0.0, scale=1.0, size=num_samples)

    # search time windows in the noise trace
    window_list = search_windows(trace=noise_trace, 
                                 threshold=num_threshold*np.std(noise_trace), 
                                 filter='off')

    return(len(window_list))

# main simulation function for all thresholds
def sim(threshold_list, num_sim, PNGname):
    transient_rate_list = []
    for num_threshold in threshold_list:
        average_windows = sum(noise_sim(num_threshold) for _ in range(num_sim)) / num_sim

        # convert average number of detected windows to transient rate
        transient_rate = average_windows / (num_samples * time_step) * 1e6 # GHz --> kHz

        print(f'\n When threshold={num_threshold},')
        print(f'\n the average transient rate over {num_sim} simulations is {transient_rate}.')
        transient_rate_list.append(transient_rate)

    # plot the relation between number of thresholds and transient rates
    plt.figure(dpi=256)
    plt.plot(threshold_list, transient_rate_list, marker='o', linestyle='-')
    plt.xlabel('Number of Threshold')
    plt.ylabel('Transient rate / kHz')
    plt.title(f'Relation between Threshold and Transient Rate \n Gaussian Noise, Separation = {standard_separation}')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'plot/{PNGname}.png')

# simulation 1
sim(np.arange(1.0, 6.0, 0.1), num_sim, 'relation1')
'''
# simulation 2
num_sim *= 2
sim(np.arange(3.14, 4.44, 0.04), num_sim, 'relation2')

# simulation 3
num_sim *= 4
sim(np.arange(3.48, 4.50, 0.01), num_sim, 'relation3')

# simulation 4
num_sim *= 5
sim(np.arange(3.79, 4.27, 0.01), num_sim, 'relation4')
'''