#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

# create the save directory if it does not exist
save_dir = os.path.join(result_dir, f'{wanted}-trace')
os.makedirs(save_dir, exist_ok=True)

# define global variables
num_abnormal = 0
num_trace    = 0

##################
# CORE FUNCTIONS #
##################

# get the traces and corresponding PSDs for 1 entry
def get1entry_trace(num_run, entry, tadc, trawv):
    global num_abnormal
    global num_trace
    be_normal = True

    # make dictionaries for easier indexing and initiate them
    filter_traces = {channel: [] for channel in channels}
    psds          = {channel: [] for channel in channels}
    filter_psds   = {channel: [] for channel in channels}

    # get info from this entry
    tadc.get_entry(entry)
    trawv.get_entry(entry)
    du            = tadc.du_id[0] # 1 entry corresponds to only 1 DU
    cst, cst_flat = get1entry_cst(trawv=trawv)
    temperature   = trawv.gps_temp[0]
    traces        = np.array(tadc.trace_ch[0])[channel_mask]
    
    for channel_id, channel in enumerate(channels):
        filter_traces[channel] = high_pass_filter(traces[channel_id])
        psds[channel]          = get_psd(traces[channel_id])
        filter_psds[channel]   = get_psd(filter_traces[channel])

        if max(filter_psds[channel][470:590]) > 1e-7:
            be_normal = False
    
    if wanted == 'normal' and not be_normal:
        return

    elif wanted == 'abnormal':
        num_trace += 1
        if be_normal:
            return
        else:
            num_abnormal += 1

    elif wanted == 'pulse':
        num_pulse = 0
        for channel in channels:
            window_list = search_windows(trace=filter_traces[channel], num_threshold=num_threshold, filter_status='off')
            num_pulse  += len(window_list)
        if num_pulse == 0:
            return
    
    # save as a NPZ file
    save_file = os.path.join(save_dir, f'{wanted}-trace_RUN{num_run}_DU{du}_{cst_flat}.npz')
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
    date_list = get_root_dates(file_list=get_root_files(date_list=False))

    for date in date_list:
        print(f'Load ROOT files on {date}.')

        for num_run in run_list:
            file_list = get_root_files(run_list=[num_run], date_list=[date])

            for file_id, file in enumerate(file_list, 1):
                # get info from this ROOT file
                tadc        = rt.TADC(file)
                trawv       = rt.TRawVoltage(file)
                num_entries = tadc.get_number_of_entries()

                print(f'{file_id:02}/{len(file_list):02}: Loop through {num_entries} entries in {file}.')
                for entry in range(num_entries): # 1 event corresponds to several entries
                    get1entry_trace(num_run=num_run, entry=entry, tadc=tadc, trawv=trawv)
    
    #
    if trace_wanted == 'abnormal':
        print(f'Number of Very Long Traces = {num_abnormal}')
        print(f'Number of All Traces = {num_trace}')

if __name__ == "__main__":
    with record_run_time():
        main()