#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

# create the save directory if it does not exist
save_dir = os.path.join(result_dir, f'{trace_wanted}-rate')
os.makedirs(save_dir, exist_ok=True)

##################
# CORE FUNCTIONS #
##################

# get transient/pulse windows for 1 day
def get1day_window(date):
    print(f'Load ROOT files on {date}.')

    for num_run in run_list:
        file_list = get_root_files(run_list=[num_run], date_list=[date])

        du_list = get_roots_dus(file_list=file_list)

        # make dictionaries for easier indexing and initiate them
        filter_windows_list = {du: {} for du in du_list}
        cst_list            = {}
        for du in du_list:
            for channel in channels:
                filter_windows_list[du][channel] = []
            cst_list[du] = []

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
                du              = tadc.du_id[0] # 1 entry corresponds to only 1 DU
                cst, cst_flat   = get1entry_cst(trawv=trawv)
                temperature     = trawv.gps_temp[0]
                traces          = np.array(tadc.trace_ch[0])[channel_mask]

                for channel_id, channel in enumerate(channels):
                    filter_trace = high_pass_filter(traces[channel_id])

                    filter_window_list = search_windows_test(trace=filter_trace, num_threshold=num_threshold, filter_status='off')
                    
                    filter_windows_list[du][channel].append(filter_window_list)
                cst_list[du].append(cst)

        # save as NPZ files
        for du in du_list:
            save_file = os.path.join(save_dir, f'{trace_wanted}-rate_RUN{num_run}_DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_cutoff{cutoff_frequency}_{date}.npz')
            np.savez(save_file,
                     **{f'filter_windows_list{channel}': filter_windows_list[du][channel] for channel in channels},
                     **{f'cst_list{channel}': cst_list[du] for channel in channels})
            print(f'Saved: {save_file}')

#################
# MAIN FUNCTION #
#################

def main():
    date_list = get_root_dates(file_list=get_root_files(date_list=False))

    for date in date_list:
        get1day_window(date=date)

if __name__ == "__main__":
    with record_run_time():
        main()