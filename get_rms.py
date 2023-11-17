#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

#########################
# GET FILES AND DU LIST #
#########################

# get selected files
file_list = sorted(glob.glob(rms_data_files))
num_files = len(file_list)

# create an empty set to store unique DUs
du_set = set()

# loop through all files to get used DUs
for file in file_list:
    tadc  = rt.TADC(file) # initiate TADC tree of this file
    tadc.get_entry(0) # get the entry from this file
    du_set.update(tadc.get_list_of_all_used_dus()) # update the set with all used DUs from this file

# convert the set to a list in order
du_list = sorted(list(du_set))

# only consider good DUs
du_list = [du for du in du_list if du in good_du_list]

# print the list of all used DUs
print(f'\nROOT files contain data from following DUs: \n{du_list}\n')

#####################################
# COMPUTE RMS DURING DIFFERENT TIME #
#####################################

# extract pure noise trace by excluding the time windows from the original trace
def extract_noise(trace, window_list):
    noise_trace = np.copy(trace)
    for start_id, stop_id in window_list:
        noise_trace[start_id:stop_id+1] = 0
        noise_trace = noise_trace[noise_trace != 0]
    return noise_trace

# make dictionary and loop through all DUs to initiate it with empty lists
noise_list = {channel: {} for channel in channels}
for du in du_list:
        for channel in channels:
            noise_list[channel][du] = []

# loop through all files
for file_id, file in enumerate(file_list[:5]):
    # get info of this file
    tadc        = rt.TADC(file)
    num_entries = tadc.get_number_of_entries()

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'Loop through the {file_id+1}/{num_files} file with {num_entries} entries: \n{file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)

        # 1 entry corresponds to only 1 DU; only consider good DUs
        du = tadc.du_id[0]
        if du not in du_list:
            continue
    
        # get traces and search for noise trace
        traces = np.array(tadc.trace_ch[0])[channel_mask]
        for channel_id, channel in enumerate(channels):
            filtered_trace = high_pass_filter(traces[channel_id])
            window_list = search_windows(filtered_trace, num_threshold*noises[channel][du], standard_separation)
            noise_trace = extract_noise(filtered_trace, window_list)
            noise_list[channel][du].extend(noise_trace)

# loop through all DUs to save
for du in du_list:
    # save all information of a DU into a NPZ file
    rms_npz_name = f'DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_cutoff{cutoff_frequency}.npz'
    rms_npz_file = os.path.join(rms_result_dir, rms_npz_name)
    np.savez(rms_npz_file, **{f'noise_list{channel}': noise_list[channel][du] for channel in channels})
    print(f'Saved: {rms_npz_file}')