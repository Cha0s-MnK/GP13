#########################
# CONFIGURE ENVIRONMENT #
#########################

# import everything from the config module
from config import *

#################
# MAIN FUNCTION #
#################

def main():
    # get JSON files
    print(f'\nLoad JSON files from RUN{num_run}.\n')
    file_list = sorted(glob(os.path.join(noise_result_dir, f'RUN{num_run}', f'*RUN{num_run}_*.json')))

    # make dictionaries for easier indexing
    noises_list = {channel: {} for channel in channels}
    du_noises   = {channel: {} for channel in channels}
    du_set      = set()

    # loop through JSON files
    for file in file_list:
        # load the data from a JSON file
        with open(file, 'r') as file:
            json_file = json.load(file)
        du_list    = json_file['du_list']
        bkg_noises = json_file['bkg_noises']
        for du in du_list:
            du_set.add(du)
            for channel in channels:
                # keys are always strings in JSON
                noises_list[channel].setdefault(du, []).append(bkg_noises[channel][str(du)])
    
    # convert the set to a sorted list
    du_list = sorted(du_set)
    for du in du_list:
        for channel in channels:
            print(f'Noises for DU{du} channel {channel}: {noises_list[channel][du]}')
            du_noises[channel][du] = np.mean(noises_list[channel][du])
    
    # compute and save the background noises to a JSON file
    noise_ref_file = os.path.join(noise_result_dir, f'noise_RUN{num_run}.json')
    with open(noise_ref_file, 'w') as json_file:
        json.dump(du_noises, json_file, indent=4)
    print(f'\nSaved: {noise_ref_file}\n')

if __name__ == "__main__":
    with record_run_time():
        main()