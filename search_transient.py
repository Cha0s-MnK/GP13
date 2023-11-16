#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

# input dates manually
date_list = ['20231027']

#############################
# SEARCH IN A SPECIFIC DATE #
#############################

def search1day(specific_date):
    print(f'\nLoad data from {data_dir}{specific_date} and search for time windows...\n')

    ##############################
    # GET ROOT FILES AND DU LIST #
    ##############################

    # update the file list according to current date
    file_list = sorted(glob.glob(os.path.join(data_dir, specific_date, '*.root')))
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

    # make dictionaries to avoid looping through all files more than once
    windows_list = {channel: {} for channel in channels}
    gps_times_list = {}

    # loop through all DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            windows_list[channel][du] = []
        gps_times_list[du] = []

    # loop through all files for the current date
    for file_id, file in enumerate(file_list, 1):
        # get info of this file
        tadc        = rt.TADC(file)
        trawv       = rt.TRawVoltage(file)
        num_entries = tadc.get_number_of_entries()

        # loop through all entries in this file (1 event corresponds to several entries)
        print(f'Loop through the {file_id}/{num_files} file with {num_entries} entries: \n{file}')
        for entry in range(num_entries):
            # get info of this entry
            tadc.get_entry(entry)
            trawv.get_entry(entry)

            # 1 entry corresponds to only 1 DU
            du = tadc.du_id[0]
            if du not in du_list:
                continue
    
            # get traces and search transients
            traces = np.array(tadc.trace_ch[0])[channel_mask]
            for channel_id, channel in enumerate(channels):
                window_list = search_windows(trace=traces[channel_id], 
                                             threshold=num_threshold*noises[channel][du], 
                                             filter=filter_state)
                windows_list[channel][du].append(window_list)
            gps_times_list[du].append(trawv.gps_time[0])

    # loop through all DUs
    for du in du_list:
        # save all information into a NPZ file
        search_result_name = f'DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_sample{num_samples}_frequency{sample_frequency}_cutoff{cutoff_frequency}_date{specific_date[:8]}.npz'
        search_result_file = os.path.join(search_result_dir, specific_date, search_result_name)
        np.savez(search_result_file,
                 # convert lists to object arrays
                 **{f'windows_list{channel}': np.array(windows_list[channel][du], dtype=object) for channel in channels},
                 gps_times=gps_times_list[du])
        print(f'Saved: {search_result_file}')

#################
# MAIN FUNCTION #
#################

def main():
    # loop through all dates
    for specific_date in date_list:
        search1day(specific_date)
    pass

if __name__ == "__main__":
    with record_run_time():
        main()