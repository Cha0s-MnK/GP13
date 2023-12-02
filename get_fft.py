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

def get_1_abnormal_psd(file):
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

    # get information of this file
    tadc        = rt.TADC(file)
    num_entries = tadc.get_number_of_entries()

    # loop through all entries in this file (1 event corresponds to several entries)
    print(f'Loop through {num_entries} entries in {file}.')
    for entry in range(num_entries):
        # get information of this entry
        tadc.get_entry(entry)

        # 1 entry corresponds to only 1 DU
        du = tadc.du_id[0]
        
        # get the traces
        traces = np.array(tadc.trace_ch[0])[channel_mask]

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

def get_pulse_fft(file):
    # get DUs of this ROOT file
    du_list = get_root_du(file)

    # get date and time from the ROOT file
    date_time, datetime_flat = get_root_datetime(file)

    # make dictionaries for easier indexing
    sum_pulse_fft    = {channel: {} for channel in channels}
    sum_noise_fft    = {channel: {} for channel in channels}
    num_used_entries = {}

    # loop through all DUs to initiate dictionaries with empty lists
    for du in du_list:
        for channel in channels:
            sum_pulse_fft[channel][du] = 0
            sum_noise_fft[channel][du] = 0
        num_used_entries[du] = 0

    # get info of this file
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
        num_used_entries[du] += 1
        
        # get the traces
        traces = np.array(tadc.trace_ch[0])[channel_mask]

        # get the GPS time and convert it to DunHuang local time and flatten it
        gps_time = trawv.gps_time[0]
        utc = gps2utc1(gps_time)
        DunHuang_time = utc + timedelta(hours=8)
        DunHuang_time_flat = DunHuang_time.strftime('%Y%m%d%H%M%S')
            
        # get FFTs, filtered FFTs and sum them up
        for id, channel in enumerate(channels):
            pulse_trace = get_pulse_trace(traces[id])
            pulse_fft_psd = get_psd(pulse_trace)
            noise_trace = get_noise_trace(traces[id])
            noise_fft_psd = get_psd(noise_trace)
            # sum up the normalized PSD for later averaging
            sum_pulse_fft[channel][du] += pulse_fft_psd
            sum_noise_fft[channel][du] += noise_fft_psd

    # loop through all DUs
    for du in du_list:
        # compute mean FFT PSDs and save all information into a NPZ file
        fft_result_file = os.path.join(fft_result_dir, f'pulseFFT_DU{du}_frequency{sample_frequency}_cutoff{cutoff_frequency}_{datetime_flat}.npz')
        np.savez(fft_result_file,
                 **{f'mean_pulse_fft{channel}': sum_pulse_fft[channel][du] / num_used_entries[du] for channel in channels},
                 **{f'mean_noise_fft{channel}': sum_noise_fft[channel][du] / num_used_entries[du] for channel in channels})
        print(f'Saved: {fft_result_file}')

    pass

def get_pulse_trace(trace):
    # apply the high-pass filter
    trace = high_pass_filter(trace=trace, 
                             sample_frequency=sample_frequency, 
                             cutoff_frequency=cutoff_frequency)

    # search for the transient/pulse windows
    window_list = search_windows(trace=trace, 
                                 threshold=num_threshold*np.std(trace))

    # initialize pulse_trace as a zero array with the same shape as trace
    pulse_trace = np.zeros_like(trace)

    # use a mask to include signal within transient/pulse windows
    mask = np.full(trace.shape, False, dtype=bool)
    for window in window_list:
        if window != []:
            start_id, stop_id = window[0], window[1]
            mask[start_id:stop_id + 1] = True

    pulse_trace[mask] = trace[mask]

    return pulse_trace

def get_noise_trace(trace):
    # apply the high-pass filter
    trace = high_pass_filter(trace=trace, 
                             sample_frequency=sample_frequency, 
                             cutoff_frequency=cutoff_frequency)

    # search for the transient/pulse windows
    window_list = search_windows(trace=trace, 
                                 threshold=num_threshold*np.std(trace))

    # initialize pulse_trace as a zero array with the same shape as trace
    noise_trace = np.zeros_like(trace)

    # use a mask to include signal within transient/pulse windows
    mask = np.full(trace.shape, True, dtype=bool)
    for window in window_list:
        if window != []:
            start_id, stop_id = window[0], window[1]
            mask[start_id:stop_id + 1] = False

    noise_trace[mask] = trace[mask]

    return noise_trace

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

# get mean FFTs of ROOT files
def main():
    # get ROOT files
    file_list = [
    'data/RUN92/GRAND.TEST-RAW-10s-ChanXYZ_20dB_12DUs_RUN92_test.20231116015403.010_dat.root',
    'data/RUN92/GRAND.TEST-RAW-10s-ChanXYZ_20dB_12DUs_RUN92_test.20231116140503.005_dat.root'
    ]

    # compute and store mean FFTs
    for file in file_list:
        get_fft(file)

if __name__ == "__main__":
    with record_run_time():
        main()