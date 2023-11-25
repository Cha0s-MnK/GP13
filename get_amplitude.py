#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

##################################
# GET AMPLITUDES FOR A ROOT FILE #
##################################

def get_file_amplitude(file):
    # get used DUs from this ROOT file
    du_list = get_root_du(file)

    # get date and time from this ROOT file
    date_time, datetime_flat = get_root_datetime(file)

    # make dictionaries for easier indexing
    amplitude_lists          = {channel: {} for channel in channels}
    filtered_amplitude_lists = {channel: {} for channel in channels}

    # loop through used DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            amplitude_lists[channel][du] = []
            filtered_amplitude_lists[channel][du] = []

    # get information of this file
    tadc        = rt.TADC(file)
    num_entries = tadc.get_number_of_entries()

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'Loop through {num_entries} entries in {file}')
    for entry in range(num_entries):
        # get information of this entry
        tadc.get_entry(entry)

        # 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
        
        # get the traces
        traces = np.array(tadc.trace_ch[0])[channel_mask]
        
        # add amplitudes to the lists
        for id, channel in enumerate(channels):
            amplitude_lists[channel][du].extend(traces[id])
            filtered_amplitude_lists[channel][du].extend((high_pass_filter(trace=traces[id], 
                                                                           sample_frequency=sample_frequency,
                                                                           cutoff_frequency=cutoff_frequency)))

    # loop through used DUs
    for du in du_list:
        # save all information into a NPZ file
        amplitude_npz_file = os.path.join(amplitude_npz_dir, f'amplitude_DU{du}_{datetime_flat}.npz')
        np.savez(amplitude_npz_file, 
                 **{f'amplitude_lists{channel}': amplitude_lists[channel][du] for channel in channels},
                 **{f'filtered_amplitude_lists{channel}': filtered_amplitude_lists[channel][du] for channel in channels})
        print(f'Saved: {amplitude_npz_file}')

    pass

#################
# MAIN FUNCTION #
#################

# get amplitudes of ROOT files
def main():
    # get the string of ROOT files
    files_str = 'data/RUN92/*.root'

    # get ROOT files
    file_list = sorted(glob(files_str))

    # get amplitudes for a file
    for file in file_list:
        get_file_amplitude(file)

if __name__ == "__main__":
    with record_run_time():
        main()