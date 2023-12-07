#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

# create the save directory if it does not exist
save_dir = os.path.join(result_dir, f'{trace_wanted}-rms-psd')
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

# get RMS and PSD from 115MHz to 144MHz of traces for 1 day
def get1day_rms_psd(date):
    print(f'Load ROOT files from {date}.')
    # loop through run numbers
    for num_run in run_list:
        # update ROOT files according to current date and run number
        file_list = get_root_files(run_list=[num_run], date=date)

        # get used DUs
        du_list = get_roots_dus(file_list=file_list)

        # make dictionaries for easier indexing
        rmses_list     = {du: {} for du in du_list}
        psd_means_list = {du: {} for du in du_list}
        psd_bands_list = {du: {} for du in du_list}
        times_list     = {du: {} for du in du_list}

        # loop through used DUs to initiate dictionaries with empty lists
        for du in du_list:
            for channel in channels:
                rmses_list[du][channel]     = []
                psd_means_list[du][channel] = []
                psd_bands_list[du][channel] = []
                times_list[du][channel]     = []

        # loop through listed files
        for file_id, file in enumerate(file_list, 1):
            # get information of this file
            tadc        = rt.TADC(file)
            trawv       = rt.TRawVoltage(file)
            num_entries = tadc.get_number_of_entries()

            # loop through all entries in this file (1 event corresponds to several entries)
            print(f'{file_id:02}/{len(file_list):02}: Loop through {num_entries} entries in {file}.')
            for entry in range(num_entries):
                # get information of this entry
                tadc.get_entry(entry)
                trawv.get_entry(entry)
                du              = tadc.du_id[0] # 1 entry corresponds to only 1 DU
                cst, cst_flat   = get1entry_cst(trawv=trawv)
                gps_temperature = trawv.gps_temp[0]
                traces          = np.array(tadc.trace_ch[0])[channel_mask]

                # loop through 3 ADC channels
                for channel_id, channel in enumerate(channels):
                    # apply the high-pass filter
                    filter_trace = high_pass_filter(traces[channel_id])

                    # get the filtered PSD
                    filter_psd = get_psd(filter_trace)

                    if trace_wanted == 'all': # store all traces
                        rmses_list[du][channel].append(rms(filter_trace))
                        psd_means_list[du][channel].append(np.mean(filter_psd[470:590]))
                        psd_bands_list[du][channel].append(filter_psd[470:590])
                        times_list[du][channel].append(cst)

        # loop through used DUs
        for du in du_list:
            # save all information into a NPZ file
            save_file = os.path.join(save_dir, f'{trace_wanted}-rms-psd_RUN{num_run}_DU{du}_{date}.npz')
            np.savez(save_file,
                     **{f'rmses_list{channel}': rmses_list[du][channel] for channel in channels},
                     **{f'psd_means_list{channel}': psd_means_list[du][channel] for channel in channels},
                     **{f'psd_bands_list{channel}': psd_bands_list[du][channel] for channel in channels},
                     **{f'times_list{channel}': times_list[du][channel] for channel in channels})
            print(f'Saved: {save_file}')

#################
# MAIN FUNCTION #
#################

def main():
    # get ROOT files
    file_list = get_root_files()

    # get the dates
    date_list = get_roots_dates(file_list=file_list)
    date_list = ['20231122', '20231123', '20231124']

    # loop through the dates
    for date in date_list:
        # get RMS and PSD from 115MHz to 144MHz of traces for 1 day
        get1day_rms_psd(date=date)

if __name__ == "__main__":
    with record_run_time():
        main()