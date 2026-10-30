import os
from pathlib import Path
import sys

base_path = next(
    p for p in [Path.cwd().resolve(), *Path.cwd().resolve().parents]
    if (p / "CODE").is_dir() and (p / "DATA").is_dir()
)
sys.path.insert(0, str(Path(base_path) / 'CODE'))
from utils.helper_functions import parcellate


if __name__ == "__main__":

    atlasfile = "/path/to/external/data/atlas/resampled_atlas_Glasser360.nii.gz" # atlas file path
    task_names = ["0back", "1back", "2back"]
    task_ids = ["01", "03", "05"]

    for subid in range(1, 36):
        for task in range(len(task_names)):
            datafile = (
                f"/path/to/external/data/residuals/FSF_files/"
                f"sub{subid:02d}_{task_names[task]}.feat/stats/res4d.nii.gz"
            )

            if not os.path.exists(datafile):
                print(f"Data file not found: {datafile}")
                continue

            fname2save = (
                f"/path/to/external/data/parcellated_residuals/"
                f"sub{subid:02d}/task_{task_ids[task]}.mat"
            )

            parcellate(atlasfile, datafile, fname2save, 360)