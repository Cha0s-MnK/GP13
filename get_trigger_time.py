from config import *


def save_reconstruction(file_list, coin_time, indx):
    with open(os.path.join(result_dir,'coincidence',f'trigger_{coin_time}.txt'), 'a') as trigger_txt,\
        open(os.path.join(result_dir,'coincidence',f'position_{coin_time}.txt'),'a') as position_txt:
        trigger_txt.write('event_id\tDU_number\ttrigger_time\n')
        position_txt.write('evet_id\tDU_number\tposition_x\tposition_y\tposition_z\n')
        du_positions_file = '/pbs/home/y/yli/data/du_position.txt'
        du_positions_dict = read_du_positions(du_positions_file)
        for file in file_list:
            du       = int(os.path.basename(file).split('_')[2][2:])
            du       = str(du)
            gps_time = np.load(file, allow_pickle = True)['gps_time']
            gps_nanosecond = np.load(file,allow_pickle = True)['gps_nanosecond']
            trace_x  = np.load(file, allow_pickle = True)['filter_tracesX']
            trace_y  = np.load(file, allow_pickle = True)['filter_tracesY']
            trace_z  = np.load(file, allow_pickle = True)['filter_tracesZ']
            window_x = search_windows(trace=trace_x, filter_status='off')
            window_y = search_windows(trace=trace_y, filter_status='off')
            window_z = search_windows(trace=trace_z, filter_status='off')
            position_x = du_positions_dict[du][0]
            position_y = du_positions_dict[du][1]
            position_z = du_positions_dict[du][2]
            if window_x == [] and window_y == [] and window_z == []:
                continue
            trigger_list = []
            for channel in ['x','y','z']:
                if eval(f'window_{channel}') != []:
                    ind_trigger = np.argmax(eval(f'trace_{channel}'))
                    trigger_list.append(ind_trigger)
            avr_time = np.average(trigger_list)*2+gps_nanosecond
            trigger_time = '{:.9f}'.format(gps_time+avr_time*1e-9)
            trigger_txt.write(f'{indx+1}\t{du}\t{trigger_time}\n')
            position_txt.write(f'{indx+1}\t{du}\t{position_x}\t{position_y}\t{position_z}\n')

    print(f'Saved : {trigger_txt}')
    print(f'Saved : {position_txt}')
    
    
    
def main():
    # get NPZ files
    print(f'\nLoad NPZ files from RUN{num_run}.\n')
    
    coin_time_list =  [20231206002842, 20231206011832, 20231206021552, 20231206022422, 20231206043422]
    for indx in range(len(coin_time_list)):
        coin_time = str(coin_time_list[indx])
        #file_list = sorted(glob(os.path.join(result_dir, f'{trace_wanted}-trace', f'{trace_wanted}-trace*{coin_time}.npz')))
        file_list = glob(os.path.join(result_dir,f'{trace_wanted}-trace',f'*{coin_time}.npz'))
        num_files = len(file_list)

        # get used DUs
        du_list = get_npz_dus(file_list=file_list)

        # loop through listed files
        print(f'\nLoop through {num_files} NPZ files from RUN{num_run}.\n')
        #for file in file_list:
        #    plot1npz_trace(file)
        save_reconstruction(file_list,coin_time,indx)
        
        
if __name__ == "__main__":
    with record_run_time():
        main()
        