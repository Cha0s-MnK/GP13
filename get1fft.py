#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

#######################################
# COMPUTE AND STORE MEAN FFT FOR A DU #
#######################################

def compute_fft(file, du_list, date_time):
    print(f'\nCompute and store FFT.\n')

    # make dictionaries for easier indexing
    sum_fft           = {channel: {} for channel in channels}
    sum_filtered_fft  = {channel: {} for channel in channels}
    mean_fft          = {channel: {} for channel in channels}
    mean_filtered_fft = {channel: {} for channel in channels}
    num_du_entries    = {}

    # loop through all DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            sum_fft[channel][du] = 0
            sum_filtered_fft[channel][du] = 0
            mean_fft[channel][du] = 0
            mean_filtered_fft[channel][du] = 0
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

    # loop through all DUs to compute mean FFTs
    for du in du_list:
        for channel in channels:
            mean_fft[channel][du]          = sum_fft[channel][du] / num_du_entries[du]
            mean_filtered_fft[channel][du] = sum_filtered_fft[channel][du] / num_du_entries[du]

        # save all information into a NPZ file
        fft_result_name = f'fft_DU{du}_threshold{num_threshold}_separation{standard_separation}_crossing{num_crossings}_max{max_samples}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{date_time}.npz'
        fft_result_file = os.path.join(fft_result_dir, fft_result_name)
        np.savez(fft_result_file,
                **{f'mean_fft{channel}': mean_fft[channel][du] for channel in channels},
                **{f'mean_filtered_fft{channel}': mean_filtered_fft[channel][du] for channel in channels})
        print(f'Saved: {fft_result_file}')

#################
# MAIN FUNCTION #
#################

def main():
    # get mean FFT of only  1 ROOT file
    fft_data_file = 'data/20231028/GRAND.TEST-RAW-10s-ChanXYZ_20dB_11DUs_RUN80_test.20231028001953.141_dat.root'

    # get date and time from the ROOT file
    date_time = get_root_datetime(fft_data_file)
    
    # get ROOT files and DUs
    file_list, du_list = get_root_du(fft_data_file)

    # get frequencies of the FFT; [MHz]
    fft_frequency = np.fft.rfftfreq(num_samples) * sample_frequency

    # compute and store mean FFTs
    compute_fft(fft_data_file, du_list, date_time)

if __name__ == "__main__":
    with record_run_time():
        main()