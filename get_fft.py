#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

#########################################
# COMPUTE AND STORE MEAN FFT FOR A FILE #
#########################################

def compute1fft(file):
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
            sum_fft[channel][du] = 0
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

        # 1 entry corresponds to only 1 DU; only consider good DUs
        du = tadc.du_id[0]
        if du not in du_list:
            continue
            
        # count this entry
        num_du_entries[du] += 1
        
        # get traces
        traces = np.array(tadc.trace_ch[0])[channel_mask]
            
        # get FFTs, filtered FFTs and sum them up
        for id, channel in enumerate(channels):
            fft = np.abs(np.fft.rfft(traces[id]))
            filtered_fft = np.abs(np.fft.rfft(high_pass_filter(trace=traces[id], 
                                                               sample_frequency=sample_frequency,
                                                               cutoff_frequency=cutoff_frequency)))
            sum_fft[channel][du] += fft
            sum_filtered_fft[channel][du] += filtered_fft

    # loop through all DUs to compute mean FFTs and save all information into a NPZ file
    for du in du_list:
        fft_result_dir  = f'result/fft/DU{du}/'
        fft_result_name = f'fft_DU{du}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{datetime_flat}.npz'
        fft_result_file = os.path.join(fft_result_dir, fft_result_name)
        np.savez(fft_result_file,
                **{f'mean_fft{channel}': sum_fft[channel][du] / num_du_entries[du] for channel in channels},
                **{f'mean_filtered_fft{channel}': sum_filtered_fft[channel][du] / num_du_entries[du] for channel in channels})
        print(f'Saved: {fft_result_file}')

    pass

#################
# MAIN FUNCTION #
#################

# get mean FFTs of only 1 ROOT file
def main():
    # get ROOT files
    file_list = [
    'data/20231014/GRAND.TEST-RAW-10s-ChanXYZ_20dB_DU10_16_17_19_20_21_29_32_33_35_10Dus_test.20231014131923.027_dat.root',
    'data/20231028/GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231028121653.161_dat.root',
    'data/20231014/GRAND.TEST-RAW-10s-ChanXYZ_20dB_DU10_16_17_19_20_21_29_32_33_35_10Dus_test.20231014015122.021_dat.root',
    'data/20231028/GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231028001953.141_dat.root'
    ]

    # compute and store mean FFTs
    for file in file_list:
        compute1fft(file)

if __name__ == "__main__":
    with record_run_time():
        main()