#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

#########################################
# SEARCH FOR TRANSIENTS/PULSES IN A DAY #
#########################################

def search_1_day(file_list):
    # get used DUs from these ROOT files
    du_list = get_root_dus(file_list)

    # get the date from these ROOT files
    date_time, datetime_flat = get_root_datetime(file_list[0])
    date = datetime_flat[:8]

    # make dictionaries for easier indexing
    windows_list          = {channel: {} for channel in channels}
    filtered_windows_list = {channel: {} for channel in channels}
    gps_time_list         = {}

    # loop through used DUs to initiate the dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            windows_list[channel][du]          = []
            filtered_windows_list[channel][du] = []
        gps_time_list[du] = []

    # loop through the files
    num_files = len(file_list)
    print(f'\nLoop through {num_files} ROOT files from {date}.\n')
    for file_id, file in enumerate(file_list, 1):
        # get information of this file
        tadc        = rt.TADC(file)
        trawv       = rt.TRawVoltage(file)
        num_entries = tadc.get_number_of_entries()

        # loop through all entries in this file (1 event corresponds to several entries)
        print(f'{file_id:02}/{num_files:02}: Loop through {num_entries} entries in {file}.')
        for entry in range(num_entries):
            # get information of this entry
            tadc.get_entry(entry)
            trawv.get_entry(entry)

            # 1 entry corresponds to only 1 DU
            du = tadc.du_id[0]
    
            # get the traces
            traces = np.array(tadc.trace_ch[0])[channel_mask]

            # search for the transient/pulses
            for channel_id, channel in enumerate(channels):
                window_list          = search_windows(trace=traces[channel_id], 
                                                      threshold=num_threshold*noises[channel][str(du)], 
                                                      filter='off')
                filtered_window_list = search_windows(trace=traces[channel_id], 
                                                      threshold=num_threshold*noises[channel][str(du)], 
                                                      filter='on')
                windows_list[channel][du].append(window_list)
                filtered_windows_list[channel][du].append(filtered_window_list)
            gps_time_list[du].append(trawv.gps_time[0])

    # create the save directory if it does not exist
    os.makedirs(os.path.join(save_dir, date), exist_ok=True)

    # loop through used DUs
    for du in du_list:
        # save all information into a NPZ file
        search_result_file = os.path.join(save_dir, date, f'search_DU{du}_RUN{num_run}_threshold{num_threshold}' \
                                                                   f'_separation{standard_separation}' \
                                                                   f'_crossing{num_crossings}_max{max_samples}' \
                                                                   f'_fluctuation{std_fluctuation}' \
                                                                   f'_frequency{sample_frequency}_cutoff{cutoff_frequency}' \
                                                                   f'_{date}.npz')
        np.savez(search_result_file,
                 # convert lists to object arrays
                 **{f'windows_list{channel}': np.array(windows_list[channel][du], dtype=object) for channel in channels},
                 **{f'filtered_windows_list{channel}': np.array(filtered_windows_list[channel][du], dtype=object) for channel in channels},
                 gps_time_list=gps_time_list[du])
        print(f'Saved: {search_result_file}')

    pass

def search_1_day_Nov(file_list):
    # get used DUs from these ROOT files
    du_list = get_root_dus(file_list)

    # get the date from these ROOT files
    date_time, datetime_flat = get_root_datetime(file_list[0])
    date = datetime_flat[:8]

    # make dictionaries for easier indexing
    filtered_windows_list = {channel: {} for channel in channels}
    gps_time_list         = {}

    # loop through used DUs to initiate the dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            filtered_windows_list[channel][du] = []
        gps_time_list[du] = []

    # loop through the files
    num_files = len(file_list)
    print(f'\nLoop through {num_files} ROOT files from {date}.\n')
    for file_id, file in enumerate(file_list, 1):
        # get information of this file
        tadc        = rt.TADC(file)
        trawv       = rt.TRawVoltage(file)
        num_entries = tadc.get_number_of_entries()

        # loop through all entries in this file (1 event corresponds to several entries)
        print(f'{file_id:02}/{num_files:02}: Loop through {num_entries} entries in {file}.')
        for entry in range(num_entries):
            # get information of this entry
            tadc.get_entry(entry)
            trawv.get_entry(entry)

            # 1 entry corresponds to only 1 DU
            du = tadc.du_id[0]
    
            # get the traces
            traces = np.array(tadc.trace_ch[0])[channel_mask]

            # search for the transient/pulses
            for channel_id, channel in enumerate(channels):
                filtered_window_list = search_windows(trace=traces[channel_id], 
                                                      threshold=num_threshold*noises[channel][str(du)], 
                                                      filter='on')
                filtered_windows_list[channel][du].append(filtered_window_list)
            gps_time_list[du].append(trawv.gps_time[0])

    # create the save directory if it does not exist
    os.makedirs(os.path.join(save_dir, date), exist_ok=True)

    # loop through used DUs
    for du in du_list:
        # save all information into a NPZ file
        search_result_file = os.path.join(save_dir, date, f'searchNov_DU{du}_RUN{num_run}_threshold{num_threshold}' \
                                                                   f'_separation{standard_separation}' \
                                                                   f'_crossing{num_crossings}_max{max_samples}' \
                                                                   f'_fluctuation{std_fluctuation}' \
                                                                   f'_frequency{sample_frequency}_cutoff{cutoff_frequency}' \
                                                                   f'_{date}.npz')
        np.savez(search_result_file,
                 # convert lists to object arrays
                 **{f'filtered_windows_list{channel}': np.array(filtered_windows_list[channel][du], dtype=object) for channel in channels},
                 gps_time_list=gps_time_list[du])
        print(f'Saved: {search_result_file}')

    pass

def search_1_day_duration(file_list):
    # get used DUs from these ROOT files
    du_list = get_root_dus(file_list)

    # get the date from these ROOT files
    date_time, datetime_flat = get_root_datetime(file_list[0])
    date = datetime_flat[:8]

    # make dictionaries for easier indexing
    durations_list = {channel: {} for channel in channels}

    # loop through used DUs to initiate the dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            durations_list[channel][du] = []

    # loop through the files
    num_files = len(file_list)
    print(f'\nLoop through {num_files} ROOT files from {date}.\n')
    for file_id, file in enumerate(file_list, 1):
        # get information of this file
        tadc        = rt.TADC(file)
        num_entries = tadc.get_number_of_entries()

        # loop through all entries in this file (1 event corresponds to several entries)
        print(f'{file_id:02}/{num_files:02}: Loop through {num_entries} entries in {file}.')
        for entry in range(num_entries):
            # get information of this entry
            tadc.get_entry(entry)

            # 1 entry corresponds to only 1 DU
            du = tadc.du_id[0]
    
            # get the traces
            traces = np.array(tadc.trace_ch[0])[channel_mask]

            # search for the transient/pulses
            for channel_id, channel in enumerate(channels):
                window_list = rough_search_windows(trace=traces[channel_id], 
                                                   threshold=num_threshold*noises[channel][str(du)])
                if window_list != []:
                    for window in window_list:
                        start_id, stop_id = window[0], window[1]
                        duration = (stop_id - start_id) * time_step
                        durations_list[channel][du].append(duration)

    # create the save directory if it does not exist
    os.makedirs(os.path.join(duration_result_dir, date), exist_ok=True)

    # loop through used DUs
    for du in du_list:
        # save all information into a NPZ file
        duration_result_file = os.path.join(duration_result_dir, date, f'duration_DU{du}_RUN{num_run}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_fluctuation{std_fluctuation}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{date}.npz')
        np.savez(duration_result_file,
                 **{f'durations_list{channel}': durations_list[channel][du] for channel in channels})
        print(f'Saved: {duration_result_file}')

    pass

#################
# MAIN FUNCTION #
#################

def main():
    # get ROOT files
    print(f'\nLoad ROOT files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(search_data_dir, '*.root')))

    # get the dates
    date_list = get_root_dates(file_list)

    # loop through all dates
    for date in date_list:
        # update ROOT files according to current date
        print(f'\nLoad ROOT files from {date}.\n')
        file_list = sorted(glob(os.path.join(search_data_dir, f'*test.{date}*.root')))

        # search for the transients/pulses
        search_1_day_Nov(file_list)
    pass

if __name__ == "__main__":
    with record_run_time():
        main()