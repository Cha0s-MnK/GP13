#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

# create the save directory if it does not exist
save_dir = os.path.join(result_dir, f'{trace_wanted}-psd')
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

def get1file_psd(file, file_id, num_files):
    # get information of this ROOT file
    tadc                = rt.TADC(file)
    trawv               = rt.TRawVoltage(file)
    num_entries         = tadc.get_number_of_entries()
    du_list             = get_root_dus(file=file)
    DHtime, DHtime_flat = get_root_DHtime(file=file)

    # make dictionaries for easier indexing
    sum_fft           = {channel: {} for channel in channels}
    sum_filtered_fft  = {channel: {} for channel in channels}
    num_du_entries    = {}

    # loop through used DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            sum_psd[channel][du]        = 0
            sum_filter_psd[channel][du] = 0
        num_du_entries[du] = 0

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'{file_id:02}/{num_files:02}: Loop through {num_entries} entries in {file}')
    for entry in range(num_entries):
        # get information of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)
        du              = tadc.du_id[0] # 1 entry corresponds to only 1 DU
        if du not in check_du_list: # only care some DUs
            return
        gps_time        = trawv.gps_time[0]
        gps_temperature = trawv.gps_temp[0]
        traces          = np.array(tadc.trace_ch[0])[channel_mask]

        # convert GPS time to CST and flatten it
        cst      = gps2utc1(gps_time) + timedelta(hours=8)
        cst_flat = cst.strftime('%Y%m%d%H%M%S')

        # get PSDs, filtered ones and sum them up
        num_du_entries[du] += 1
        for channel_id, channel in enumerate(channels):
            psd        = get_psd(trace=traces[channel_id])
            filter_psd = get_psd(trace=high_pass_filter(trace=traces[channel_id]))

            # sum up the PSDs for later averaging
            sum_psd[channel][du]        += psd
            sum_filter_psd[channel][du] += filter_psd
        
    # loop through used DUs
    for du in du_list:
        # save all information into a NPZ file
        fft_result_file = os.path.join(save_dir, f'FFT_DU{du}_RUN{num_run}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{datetime_flat}.npz')
        np.savez(fft_result_file,
                **{f'mean_fft{channel}': sum_fft[channel][du] / num_du_entries[du] for channel in channels},
                **{f'mean_filtered_fft{channel}': sum_filtered_fft[channel][du] / num_du_entries[du] for channel in channels})
        print(f'Saved: {fft_result_file}')

def get_1_abnormal_psd(file):
        for channel_id, channel in enumerate(channels):
            # check if this trace has an abnormal fluctuation of standard deviation
            if np.std(traces[channel_id]) < (std_fluctuation * noises[channel][str(du)]):
                # get FFT PSD and filtered FFT PSD
                for id, channel in enumerate(channels):
                    fft_psd          = get_psd(traces[id])
                    filtered_fft_psd = get_psd((high_pass_filter(trace=traces[id], 
                                                                 sample_frequency=sample_frequency,
                                                                 cutoff_frequency=cutoff_frequency)))

    # loop through all DUs
    for du in du_list:
        # create the saved directory if it does not exist
        os.makedirs(os.path.join(fft_result_dir, f'DU{du}'), exist_ok=True)

        # compute mean FFT PSDs and save all information into a NPZ file
        fft_result_file = os.path.join(fft_result_dir, f'FFT_DU{du}_RUN{num_run}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{datetime_flat}.npz')
        np.savez(fft_result_file,
                **{f'mean_fft{channel}': sum_fft[channel][du] / num_du_entries[du] for channel in channels},
                **{f'mean_filtered_fft{channel}': sum_filtered_fft[channel][du] / num_du_entries[du] for channel in channels})
        print(f'Saved: {fft_result_file}')

    pass


def get_fft(file):
    # get DUs of this ROOT file
    du_list = get_root_du(file)

    # get date and time from the ROOT file
    date_time, datetime_flat = get_root_datetime(file)

    # make dictionaries for easier indexing
    sum_fft           = {channel: {} for channel in channels}
    sum_filtered_fft  = {channel: {} for channel in channels}
    num_du_entries    = {}

    # loop through all DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            sum_fft[channel][du]          = 0
            sum_filtered_fft[channel][du] = 0
        num_du_entries[du] = 0

    # get info of this file
    tadc        = rt.TADC(file)
    num_entries = tadc.get_number_of_entries()

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'Loop through the file with {num_entries} entries: \n{file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)

        # 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
        num_du_entries[du] += 1
        
        # get traces
        traces = np.array(tadc.trace_ch[0])[channel_mask]
        
        # get FFTs, filtered FFTs and sum them up
        for id, channel in enumerate(channels):
            fft_psd          = get_psd(traces[id])
            filtered_fft_psd = get_psd((high_pass_filter(trace=traces[id], 
                                                         sample_frequency=sample_frequency,
                                                         cutoff_frequency=cutoff_frequency)))

            # sum up the normalized PSD for later averaging
            sum_fft[channel][du] += fft_psd
            sum_filtered_fft[channel][du] += filtered_fft_psd

    # loop through all DUs
    for du in du_list:
        # create the saved directory if it does not exist
        os.makedirs(os.path.join(fft_result_dir, f'DU{du}'), exist_ok=True)

        # compute mean FFT PSDs and save all information into a NPZ file
        fft_result_file = os.path.join(fft_result_dir, f'FFT_DU{du}_RUN{num_run}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{datetime_flat}.npz')
        np.savez(fft_result_file,
                **{f'mean_fft{channel}': sum_fft[channel][du] / num_du_entries[du] for channel in channels},
                **{f'mean_filtered_fft{channel}': sum_filtered_fft[channel][du] / num_du_entries[du] for channel in channels})
        print(f'Saved: {fft_result_file}')

    pass

#################
# MAIN FUNCTION #
#################

def main():
    # get ROOT files
    file_list = sorted(glob(os.path.join(data_dir, f'RUN{num_run}', '*.root')))
    num_files = len(file_list)
    print(f'\nLoad {num_files} ROOT files from RUN{num_run}.\n')

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
            # get mean PSDs for 1 file
            get1file_psd(file=file, file_id=file_id, num_files=num_files)

if __name__ == "__main__":
    with record_run_time():
        main()