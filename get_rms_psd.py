#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

# create the save directory if it does not exist
name = 'rms-psd'
save_dir = os.path.join(result_dir, name)
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

# get RMS and mean PSDs for 1 day
def get1day_rms_psd(date):
    print(f'Load ROOT files on {date}.')

    for num_run in run_list:
        file_list = get_root_files(run_list=[num_run], date_list=[date])
        du_list   = get_roots_dus(file_list=file_list)

        # make dictionaries for easier indexing and initiate them
        rmses_list        = {du: {} for du in du_list}
        psds_list         = {du: {} for du in du_list}
        times_list        = {du: {} for du in du_list}
        temperatures_list = {du: {} for du in du_list}
        for du in du_list:
            for channel in channels:
                rmses_list[du][channel]        = []
                psds_list[du][channel]         = []
                times_list[du][channel]        = []
                temperatures_list[du][channel] = []

        for file_id, file in enumerate(file_list, 1):
            # get info from this ROOT file
            tadc        = rt.TADC(file)
            trawv       = rt.TRawVoltage(file)
            num_entries = tadc.get_number_of_entries()

            print(f'{file_id:02}/{len(file_list):02}: Loop through {num_entries} entries in {file}.')
            for entry in range(num_entries): # 1 event corresponds to several entries
                # get info from this entry
                tadc.get_entry(entry)
                trawv.get_entry(entry)
                du            = tadc.du_id[0] # 1 entry corresponds to only 1 DU
                cst, cst_flat = get1entry_cst(trawv=trawv)
                temperature   = trawv.gps_temp[0]
                traces        = np.array(tadc.trace_ch[0])[channel_mask]

                for channel_id, channel in enumerate(channels):
                    filter_trace = high_pass_filter(traces[channel_id])
                    filter_psd = get_psd(filter_trace)

                    rmses_list[du][channel].append(rms(filter_trace))
                    psds_list[du][channel].append(filter_psd[205:820]) # 50-200 MHz
                    times_list[du][channel].append(cst)
                    temperatures_list[du][channel].append(temperature)

        # save as a NPZ file
        for du in du_list:
            save_file = os.path.join(save_dir, f'{name}_RUN{num_run}_DU{du}_{date}.npz')
            np.savez(save_file,
                     **{f'rmses_list{channel}': rmses_list[du][channel] for channel in channels},
                     **{f'psds_list{channel}': psds_list[du][channel] for channel in channels},
                     **{f'times_list{channel}': times_list[du][channel] for channel in channels},
                     **{f'temperatures_list{channel}': temperatures_list[du][channel] for channel in channels})
            print(f'Saved: {save_file}')

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

#################
# MAIN FUNCTION #
#################

def main():
    date_list = get_root_dates(file_list=get_root_files(date_list=['20231207']))

    for date in date_list:
        get1day_rms_psd(date=date)

if __name__ == "__main__":
    with record_run_time():
        main()