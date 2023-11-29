#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

######################################
# COMPUTE BACKGROUND NOISES IN A DAY #
######################################

def compute_noise(file_list):
    # get used DUs from these ROOT files
    du_list = get_root_dus(file_list)

    # get the date from these ROOT files
    date_time, datetime_flat = get_root_datetime(file_list[0])
    date = datetime_flat[:8]

    # make dictionaries for easier indexing
    noises_list = {channel: {} for channel in channels}
    bkg_noises  = {channel: {} for channel in channels}

    # loop through used DUs to initiate the dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            noises_list[channel][du] = []
            bkg_noises[channel][du]  = 0

    # loop through the files
    num_files = len(file_list)
    print(f'\nLoop through {num_files} ROOT files from {date}.\n')
    for file_id, file in enumerate(file_list, 1):
        # get information of this file
        tadc        = rt.TADC(file)
        num_entries = tadc.get_number_of_entries()

        # loop through all entries in this file (1 event corresponds to several entries)
        print(f'{file_id:02}/{num_files:02}: Loop through {num_entries} entries in {file}')
        for entry in range(num_entries):
            # get information of this entry
            tadc.get_entry(entry)

            # 1 entry corresponds to only 1 DU
            du = tadc.du_id[0]
    
            # get the traces
            traces = np.array(tadc.trace_ch[0])[channel_mask]
            for channel_id, channel in enumerate(channels):
                # apply the high-pass filter
                filtered_trace = high_pass_filter(trace=traces[channel_id], 
                                                  sample_frequency=sample_frequency, 
                                                  cutoff_frequency=cutoff_frequency)

                # primitive search for the transient/pulse windows
                window_list = search_windows(trace=filtered_trace, 
                                             threshold=num_threshold*np.std(filtered_trace))

                # exclude signals within transient/pulse windows
                mask = np.full(filtered_trace.shape, True, dtype=bool)
                for window in window_list:
                    if window != []:
                        start_id, stop_id = window[0], window[1]
                        mask[start_id:stop_id + 1] = False
                noise_trace = filtered_trace[mask]

                # consider the standard deviation of the noise trace as the background noise
                noises_list[channel][du].append(np.std(noise_trace))

    for du in du_list:
        for channel in channels:
            # exclude abnormally large background noises
            noise_array = np.array(noises_list[channel][du])
            mean = np.mean(noise_array)
            mask = noise_array < std_fluctuation * mean
            masked_noise_array = noise_array[mask]

            # 
            masked_mean = np.mean(masked_noise_array)
            bkg_noises[channel][du] = masked_mean

    # save the DUs and the background noise dictionary to a JSON file
    json_dict = {
        'du_list': du_list,
        'bkg_noises': bkg_noises
    }
    with open(os.path.join(noise_result_dir, f'noise_RUN{num_run}_{date}.json'), 'w') as json_file:
        json.dump(json_dict, json_file, indent=4)
    
    pass

#################
# MAIN FUNCTION #
#################

def main():
    # get ROOT files
    print(f'\nLoad ROOT files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(noise_data_dir, '*.root')))

    # get the dates
    date_list = get_root_dates(file_list)

    # loop through all dates
    for date in date_list:
        # update ROOT files according to current date
        print(f'\nLoad ROOT files from {date}.\n')
        file_list = sorted(glob(os.path.join(noise_data_dir, f'*test.{date}*.root')))

        # compute background noise for 1 day
        compute_noise(file_list)
    pass

if __name__ == "__main__":
    with record_run_time():
        main()