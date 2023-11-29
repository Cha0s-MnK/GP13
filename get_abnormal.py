#########################
# CONFIGURE ENVIRONMENT #
#########################

# enable GRANDLIB
import grand.dataio.root_trees as rt

# import everything from the config module
from config import *

########################################
# COMPUTE MEAN FFT PSD FOR A ROOT FILE #
########################################

def get_root_abnormal(file):
    # get information of this file
    tadc        = rt.TADC(file)
    trawv       = rt.TRawVoltage(file)
    num_entries = tadc.get_number_of_entries()

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'Loop through {num_entries} entries in {file}.')
    for entry in range(num_entries):
        # get information of this entry
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
        
        # get the traces
        traces = np.array(tadc.trace_ch[0])[channel_mask]

        # get the GPS time and convert it to DunHuang time (DHtime) and flatten it
        gps_time    = trawv.gps_time[0]
        DHtime      = gps2utc1(gps_time) + timedelta(hours=8)
        DHtime_flat = DHtime.strftime('%Y%m%d%H%M%S')

        for channel_id, channel in enumerate(channels):
            # check if this trace has an abnormal fluctuation of standard deviation
            if np.std(traces[channel_id]) > (std_fluctuation * noises[channel][str(du)]):
                # get FFT PSD and filtered FFT PSD
                for id, channel in enumerate(channels):
                    psd        = get_psd(traces[id])
                    filter_psd = get_psd((high_pass_filter(trace=traces[id], 
                                                           sample_frequency=sample_frequency,
                                                           cutoff_frequency=cutoff_frequency)))

                # save the traces into a NPZ file
                abnormal_result_file = os.path.join(abnormal_result_dir, f'abnormal_DU{check_du}_RUN{num_run}_fluctuation{std_fluctuation}_{DHtime_flat}.npz')
                np.savez(abnormal_result_file,
                         trace=traces[channel_id],
                         channel=channel,
                         psd=psd,
                         filter_psd=filter_psd)
                print(f'Saved: {abnormal_result_file}')

    pass

def get_psd(trace):
    # get FFTs
    fft = np.abs(np.fft.rfft(trace))

    # convert FFT from ADC units to Volts and adjust for system linear gain
    fft = fft * adcu2v / linear_gain

    # compute power of FFT and normalize it to number of samples in the trace
    num_samples = len(trace)
    fft_power   = fft * fft / num_samples / num_samples

    # calculate the frequency bin width
    fft_frequency       = np.fft.rfftfreq(num_samples) * sample_frequency # frequencies of the FFT [MHz]
    frequency_bin_width = fft_frequency[1] - fft_frequency[0]

    # normalize power to frequency bin width to get PSD
    fft_psd = fft_power / frequency_bin_width
    return fft_psd

#################
# MAIN FUNCTION #
#################

def main():
    # get ROOT files
    print(f'\nLoad ROOT files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(abnormal_data_dir, '*.root')))

    # loop through all files
    for file in file_list:
        get_root_abnormal(file)

if __name__ == "__main__":
    with record_run_time():
        main()