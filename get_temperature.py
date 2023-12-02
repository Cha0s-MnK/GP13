#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

##################
# CORE FUNCTIONS #
##################

# get background noises and GPS temperatures for 1 day
def get1day_temperature(date):
    # update ROOT files according to current date
    file_list = sorted(glob(os.path.join(data_dir, f'RUN{num_run}', f'*test.{date}*.root')))

    # get used DUs
    du_list = get_root_dus(file_list=file_list)

    # make dictionaries for easier indexing
    noises_list      = {channel: {} for channel in channels}
    temperature_list = {}

    # loop through used DUs to initiate the dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            noises_list[channel][du] = []
            temperature_list[du]     = []

    # loop through listed files
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
            du              = tadc.du_id[0] # 1 entry corresponds to only 1 DU
            gps_time        = trawv.gps_time[0]
            gps_temperature = trawv.gps_temp[0]
            traces          = np.array(tadc.trace_ch[0])[channel_mask]

            # convert GPS time to DunHuang time and flatten it
            DHtime      = gps2utc1(gps_time) + timedelta(hours=8)
            DHtime_flat = DHtime.strftime('%Y%m%d%H%M%S')

            # compute background noises
            for channel_id, channel in enumerate(channels):
                # apply the high-pass filter
                filter_trace = high_pass_filter(trace=traces[channel_id])

                # primitive search for transient/pulse windows
                window_list = search_windows(trace=filter_trace, threshold=num_threshold*np.std(filter_trace))

                # exclude signals within transient/pulse windows
                mask = np.full(filter_trace.shape, True, dtype=bool)
                for window in window_list:
                    if window != []:
                        start_id, stop_id = window[0], window[1]
                        mask[start_id:stop_id + 1] = False
                noise_trace = filter_trace[mask]

                # take the standard deviation of the noise trace as the background noise
                noises_list[channel][du].append(np.std(noise_trace))
            
            # add GPS temperature
            temperature_list[du].append(gps_temperature)

    # create the save directory if it does not exist
    save_dir = os.path.join(result_dir, 'temperature')
    os.makedirs(save_dir, exist_ok=True)

    # loop through used DUs
    for du in du_list:
        # save all information into a NPZ file
        save_file = os.path.join(save_dir, f'temperature_RUN{num_run}_DU{du}_{date}.npz')
        np.savez(save_file,
                 **{f'noises_list{channel}': noises_list[channel][du] for channel in channels},
                 temperature_list=temperature_list[du])
        print(f'Saved: {save_file}')

# get normal background noises and GPS temperatures for 1 day
def get1day_normal_temperature(date):
    # update ROOT files according to current date
    file_list = sorted(glob(os.path.join(data_dir, f'RUN{num_run}', f'*test.{date}*.root')))

    # get used DUs
    du_list = get_root_dus(file_list=file_list)

    # make dictionaries for easier indexing
    noises_list      = {channel: {} for channel in channels}
    temperature_list = {}

    # loop through used DUs to initiate the dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            noises_list[channel][du] = []
            temperature_list[du]     = []

    # loop through listed files
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
            du              = tadc.du_id[0] # 1 entry corresponds to only 1 DU
            gps_time        = trawv.gps_time[0]
            gps_temperature = trawv.gps_temp[0]
            traces          = np.array(tadc.trace_ch[0])[channel_mask]

            # convert GPS time to DunHuang time and flatten it
            DHtime      = gps2utc1(gps_time) + timedelta(hours=8)
            DHtime_flat = DHtime.strftime('%Y%m%d%H%M%S')

            # loop through 3 ADC channels to check whether this entry is abnormal/noisy
            be_abnormal = False
            for channel_id, channel in enumerate(channels):
                # apply the high-pass filter
                filter_trace = high_pass_filter(trace=traces[channel_id])

                # get PSD of the filtered trace
                filter_psd = get_psd(filter_trace)

                # check whether this trace is abnormal/noisy
                if max(filter_psd) > 1e-7:
                    be_abnormal = True
                    break
            
            if be_abnormal:
                print('abnormal!')
                continue

            # compute background noises
            for channel_id, channel in enumerate(channels):
                # apply the high-pass filter
                filter_trace = high_pass_filter(trace=traces[channel_id])

                # primitive search for transient/pulse windows
                window_list = search_windows(trace=filter_trace, threshold=num_threshold*np.std(filter_trace))

                # exclude signals within transient/pulse windows
                mask = np.full(filter_trace.shape, True, dtype=bool)
                for window in window_list:
                    if window != []:
                        start_id, stop_id = window[0], window[1]
                        mask[start_id:stop_id + 1] = False
                noise_trace = filter_trace[mask]

                # take the standard deviation of the noise trace as the background noise
                noises_list[channel][du].append(np.std(noise_trace))
            
            # add GPS temperature
            temperature_list[du].append(gps_temperature)

    # create the save directory if it does not exist
    save_dir = os.path.join(result_dir, 'temperature')
    os.makedirs(save_dir, exist_ok=True)

    # loop through used DUs
    for du in du_list:
        # save all information into a NPZ file
        save_file = os.path.join(save_dir, f'normal-temperature_RUN{num_run}_DU{du}_{date}.npz')
        np.savez(save_file,
                 **{f'noises_list{channel}': noises_list[channel][du] for channel in channels},
                 temperature_list=temperature_list[du])
        print(f'Saved: {save_file}')

#################
# MAIN FUNCTION #
#################

def main():
    # get ROOT files
    print(f'\nLoad ROOT files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(data_dir, f'RUN{num_run}', '*.root')))

    # get the dates
    date_list = get_root_dates(file_list=file_list)

    # get used DUs
    du_list = get_root_dus(file_list=file_list)

    # loop through the dates
    for date in date_list:
        # get background noises and GPS temperatures for 1 day
        get1day_normal_temperature(date)

if __name__ == "__main__":
    with record_run_time():
        main()