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

# get the traces and corresponding PSDs for 1 entry
def get1entry_trace(entry, tadc, trawv):
    # make dictionaries for easier indexing and initiate them with empty lists
    filter_traces = {channel: [] for channel in channels}
    psds          = {channel: [] for channel in channels}
    filter_psds   = {channel: [] for channel in channels}

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

    # loop through 3 ADC channels
    for channel_id, channel in enumerate(channels):
        # apply the high-pass filter
        filter_traces[channel] = high_pass_filter(traces[channel_id])

        # get the PSD and the filtered one
        psds[channel]        = get_psd(traces[channel_id])
        filter_psds[channel] = get_psd(filter_traces[channel])

    # create the save directory if it does not exist
    save_dir = os.path.join(result_dir, 'trace')
    os.makedirs(save_dir, exist_ok=True)
    
    # save all information into a NPZ file
    save_file = os.path.join(save_dir, f'trace_RUN{num_run}_DU{du}_{DHtime_flat}.npz')
    np.savez(save_file, 
             **{f'traces{channel}': traces[channel_id] for channel_id, channel in enumerate(channels)},
             **{f'filter_traces{channel}': filter_traces[channel] for channel in channels},
             **{f'psds{channel}': psds[channel] for channel in channels},
             **{f'filter_psds{channel}': filter_psds[channel] for channel in channels})
    print(f'Saved: {save_file}')

# get the traces and corresponding PSDs for 1 abnormal/noisy entry
def get1entry_abnormal_trace(entry, tadc, trawv):
    # make dictionaries for easier indexing and initiate them with empty lists
    filter_traces = {channel: [] for channel in channels}
    psds          = {channel: [] for channel in channels}
    filter_psds   = {channel: [] for channel in channels}

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

    # loop through 3 ADC channels
    be_abnormal = False
    for channel_id, channel in enumerate(channels):
        # apply the high-pass filter
        filter_traces[channel] = high_pass_filter(traces[channel_id])

        # get the PSD and the filtered one
        psds[channel]        = get_psd(traces[channel_id])
        filter_psds[channel] = get_psd(filter_traces[channel])

        # check filtered PSD to judge whether the trace is abnormal/noisy
        if max(filter_psds[channel]) > 1e-7:
            be_abnormal = True

    # return early if this is not an abnormal/noisy entry
    if not be_abnormal:
        return

    # create the save directory if it does not exist
    save_dir = os.path.join(result_dir, 'trace')
    os.makedirs(save_dir, exist_ok=True)
    
    # save all information into a NPZ file
    save_file = os.path.join(save_dir, f'abnormal-trace_RUN{num_run}_DU{du}_{DHtime_flat}.npz')
    np.savez(save_file, 
             **{f'traces{channel}': traces[channel_id] for channel_id, channel in enumerate(channels)},
             **{f'filter_traces{channel}': filter_traces[channel] for channel in channels},
             **{f'psds{channel}': psds[channel] for channel in channels},
             **{f'filter_psds{channel}': filter_psds[channel] for channel in channels})
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
        # update ROOT files according to current date
        file_list = sorted(glob(os.path.join(data_dir, f'RUN{num_run}', f'*test.{date}*.root')))

        # loop through listed files
        num_files = len(file_list)
        print(f'\nLoop through {num_files} ROOT files from {date}.\n')
        for file_id, file in enumerate(file_list, 1):
            # get information of this file
            tadc        = rt.TADC(file)
            trawv       = rt.TRawVoltage(file)
            num_entries = tadc.get_number_of_entries()

            # loop through all entries in this file (1 event corresponds to several entries)
            print(f'{file_id:02}/{num_files:02}: Loop through {num_entries} entries in {file}')
            for entry in range(num_entries):
                # get the traces and corresponding PSDs for 1 entry
                get1entry_trace(entry, tadc, trawv)

if __name__ == "__main__":
    with record_run_time():
        main()