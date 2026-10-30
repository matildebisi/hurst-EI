

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

    # -----------------------------------------------------------------
    # Paths
    # -----------------------------------------------------------------

    atlasfile = (
        "/path/to/atlas/atlas_Glasser360.nii.gz"
    )

    base_dir = "/path/to/WAND/data"

    session = "ses-03"

    output_dir = (
        "/path/to/"
        "WAND_parcellated/Glasser360"
    )

    data_filename = "s6w-noBPF-residuals.nii"

    expected_nreg = 100

    # IMPORTANT:
    # True = overwrite the old parcellated files
    overwrite = False

    # -----------------------------------------------------------------
    # Subjects
    # -----------------------------------------------------------------

    subjects = sorted([
        d
        for d in os.listdir(base_dir)
        if d.startswith("sub-")
        and os.path.isdir(os.path.join(base_dir, d))
    ])

    print(f"Found {len(subjects)} subjects.")

    # -----------------------------------------------------------------
    # Loop subjects
    # -----------------------------------------------------------------

    for subj in subjects:

        func_dir = os.path.join(
            base_dir,
            subj,
            session,
            "func"
        )

        datafile = os.path.join(
            func_dir,
            data_filename
        )

        # Try .nii.gz
        if not os.path.exists(datafile):

            datafile_gz = datafile + ".gz"

            if os.path.exists(datafile_gz):
                datafile = datafile_gz

            else:
                print(
                    f"[SKIP] Missing preprocessed file "
                    f"for {subj}: {datafile}"
                )
                continue

        fname2save = os.path.join(
            output_dir,
            subj,
            "s6w-noBPF-residuals.mat"
        )

        # -------------------------------------------------------------
        # Existing output
        # -------------------------------------------------------------

        if os.path.exists(fname2save) and not overwrite:

            print(
                f"[SKIP] Already parcellated: {subj}"
            )

            continue

        print("\n-----------------------------------------")
        print(f"[RUN] Parcellating {subj}")
        print(f"      data: {datafile}")
        print(f"      save: {fname2save}")

        try:

            parcellate(
                atlasfile=atlasfile,
                datafile=datafile,
                fname2save=fname2save,
                mean_center=True,
                expected_nreg=expected_nreg
            )

        except Exception as e:

            print(
                f"[ERROR] Failed for {subj}: {e}"
            )

            continue

    print("\nDone.")