from eyelinkio import read_edf
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.interpolate import PchipInterpolator
import pandas as pd
import os

#########################################
## preprocessing functions ##

def remove_blinks(blinks, data, times, margin=0.2):
    """Remove blinks from data by interpolation. """

    data_clean = data.copy()
    for _, start, end in blinks:
        blink_mask = (times >= start - margin) & (times <= end + margin)
        before = np.where(times < start - margin)[0] 
        after = np.where(times > end + margin)[0]    
        
        if before.size > 0 and after.size > 0:
            val_before = data_clean[before[-1]] 
            val_after = data_clean[after[0]]   
            t0 = times[before[-1]]
            t1 = times[after[0]]
            interp_vals = PchipInterpolator([t0, t1], [val_before, val_after])(times[blink_mask])
            data_clean[blink_mask] = interp_vals

        elif before.size > 0:
            data_clean[blink_mask] = data_clean[before[-1]]
        elif after.size > 0:
            data_clean[blink_mask] = data_clean[after[0]]
        else:
            data_clean[blink_mask] = np.nan
    
    return data_clean, times

# as in Behav Res (2019) 51:1336–1342
def speed_outliers_detection(data, times, n):

    speed = np.abs(np.diff(data) / np.diff(times))
    speed_back = np.insert(speed, 0, 0)
    speed_forward = np.insert(speed, speed.size, 0)
    d_prime = np.maximum(speed_back, speed_forward)

    MAD = np.median(np.abs(d_prime - np.median(d_prime)))
    threshold = np.median(d_prime) + n * MAD

    mask = d_prime <= threshold
    filtered_data = data.copy()

    if not np.all(mask): 
        good_times = times[mask]
        good_data = data[mask]
        interp_vals = PchipInterpolator(good_times, good_data)(times)
        filtered_data[~mask] = interp_vals[~mask]

    return filtered_data, d_prime

def butter_lowpass_filter(data, cutoff, fs, order=4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y
    
def preprocess_pupil(edf_file, margin, n):

    data, times = remove_blinks(edf_file['discrete']['blinks'], edf_file['samples'][2], edf_file['times'], margin)

    if np.nansum(data) == 0:
        print('All data is NaN or 0 after blink removal.')
        return None, None
    
    data, _ = speed_outliers_detection(data, times, n)
    data = butter_lowpass_filter(data, cutoff=9, fs=500, order=4)
    return data, times

#########################################



##########################################

if __name__ == "__main__":

    num_blocks = 4
    trial_per_block = 44
    num_trials = trial_per_block * num_blocks

    sub_ids = list(range(1, 37))
    task_ids = [1, 2, 3, 4, 5, 6]

    ### Paths ###

    data_path = "/data/eye_tracking/raw_data"
    output_path = "/data/eye_tracking/preprocessed_data/pupil_size"

    total_data = []

    for sub_id in sub_ids:

        for i, task_id in enumerate(task_ids):
            data_dict = {
            'subject':sub_id,                                       # Subject ID
            'task_ids': task_id,                                    # Task ID
            'preprocessed_timeseries': None,                        # Preprocessed pupil timeseries
            'times': None,                                          # Corresponding time points
            }
            
            file_name = f's{sub_id}00{task_id}.edf'
            edf_path = f'{data_path}/{file_name}'

            if not os.path.exists(edf_path):
                print(f'File not found: {edf_path}')
                continue

            edf_file = read_edf(edf_path)
            
            series, times = preprocess_pupil(edf_file, margin=0.1, n=2)

            if series is None:
                print(f'No data for subject {sub_id}, task {task_id}')
                continue

            data_dict['preprocessed_timeseries'] = series
            data_dict['times'] = times

            print("Processed subject", sub_id, "task", task_id)

            total_data.append(data_dict)

    df = pd.DataFrame(total_data)
    df.to_pickle(output_path)