#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

###############################
# SAVE GOOD TRACES OF GOOD DU #
###############################

def check1root(file):
    # get info of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'Loop through the file with {num_entries} entries: \n{file}')
    for entry in range(num_entries):
        # get info of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # 1 entry corresponds to only 1 DU; only consider the good DU
        if tadc.du_id[0] != good_du:
            continue

        # get the traces and judge whether there are time windows in channel X and Y
        traces = np.array(tadc.trace_ch[0])[channel_mask]
        thresholdX = num_threshold * noises['X'][good_du]
        thresholdY = num_threshold * noises['Y'][good_du]
        window_listX = search_windows(traces[0], thresholdX, 'on')
        window_listY = search_windows(traces[1], thresholdY, 'on')
        if len(window_listX) == 0 and len(window_listY) == 0:
            continue

        # get GPS time
        gps_time = trawv.gps_time[0]

        # convert GPS time to UTC and flatten it
        utc = gps2utc1(gps_time)
        utc_flat = utc.strftime('%Y%m%d%H%M%S')
        
        # save the traces into a NPZ file
        check_result_name = f'check_DU{good_du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{utc_flat}.npz'
        check_result_file = os.path.join(check_result_dir, check_result_name)
        np.savez(check_result_file,
                 **{f'trace{channel}': traces[id] for id, channel in enumerate(channels)},
                 utc=utc_flat)
        print(f'Saved: {check_result_file}')
    pass

#################
# MAIN FUNCTION #
#################

def main():
    # get ROOT files
    file_list = sorted(glob(check_data_files))

    # save the good traces of the good DU in these ROOT files
    for file in file_list:
        check1root(file)
    pass

if __name__ == "__main__":
    with record_run_time():
        main()