from eyelinkio import read_edf
import pandas as pd 
import os 

if __name__ == "__main__":

    sub_ids = list(range(1, 37))
    task_ids = [1, 2, 3, 4, 5, 6]

    data_path = "/data/eye_tracking/raw_data"
    output_path = "/data/eye_tracking/preprocessed_data/microsaccades"

    total_data = []

    for sub_id in sub_ids:
        sub_name = f'sub{sub_id:02d}'

        for i, task_id in enumerate(task_ids):
            data_dict = {
            'subject':sub_id,                                       # Subject ID
            'task_ids': task_id,                                    # Task ID
            'start_times': None,                                
            'end_times': None,                                  
            'start_x_positions': None,                             
            'start_y_positions': None,
            }
            
            file_name = f's{sub_id}00{task_id}.edf'
            edf_path = f'{data_path}/{file_name}'

            if not os.path.exists(edf_path):
                print(f'File not found: {edf_path}')
                continue

            edf_file = read_edf(edf_path)
            
            saccades = edf_file['discrete']['saccades']

            if saccades is None:
                print(f'No data for subject {sub_id}, task {task_id}')
                continue

            data_dict['start_times'] = [i[1] for i in saccades]
            data_dict['end_times'] = [i[2] for i in saccades]
            data_dict['start_x_positions'] = [i[3] for i in saccades]
            data_dict['start_y_positions'] = [i[4] for i in saccades]

            print("Processed subject", sub_id, "task", task_id)

            total_data.append(data_dict)

    df = pd.DataFrame(total_data)
    df.to_pickle(output_path)