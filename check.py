#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

#####################################
# SAVE TRACES OF ONE DU OF ONE FILE #
#####################################

def check(file):
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

        # get GPS time
        gps_time = trawv.gps_time[0]

        # convert GPS time to UTC and flatten it
        utc = gps2utc1(gps_time)
        utc_flat = utc.strftime('%Y%m%d%H%M%S')
        
        # get the traces and save it into a NPZ file
        traces = np.array(tadc.trace_ch[0])[channel_mask]
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

# check 1 ROOT file to test the search method
def main():
    # save the traces of the good DU in this file
    check('data/20231028/GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231028001953.141_dat.root')
    pass

if __name__ == "__main__":
    with record_run_time():
        main()