"""Build parcel-by-subject/task matrices and persist SVD results.

The imputation, ordering, centering and output filenames intentionally match
the replication notebook so the saved matrices and downstream figures remain
compatible with earlier analysis runs.
"""

import os
import numpy as np
import pandas as pd


def prepare_pivot_table(data_df, tasks):
    df = data_df.copy()

    df['Subject'] = pd.to_numeric(df['Subject'], errors='coerce')
    df['Parcel'] = pd.to_numeric(df['Parcel'], errors='coerce')
    df['Hursts'] = pd.to_numeric(df['Hursts'], errors='coerce')

    df = df.dropna(subset=['Subject', 'Parcel', 'Task', 'Hursts'])

    df['Subject'] = df['Subject'].astype(int)
    df['Parcel'] = df['Parcel'].astype(int)

    parcels = sorted(df['Parcel'].unique())
    subjects = sorted(df['Subject'].unique())

    df_pivot = (
        df.pivot_table(
            index=['Subject', 'Parcel'],
            columns='Task',
            values='Hursts',
            aggfunc='mean'
        )
        .reset_index()
    )

    full_index = pd.MultiIndex.from_product(
        [subjects, parcels],
        names=['Subject', 'Parcel']
    )

    df_pivot = (
        df_pivot
        .set_index(['Subject', 'Parcel'])
        .reindex(full_index)
        .reset_index()
    )

    for task in tasks:
        if task not in df_pivot.columns:
            df_pivot[task] = np.nan

    return df_pivot, subjects, parcels


def impute_missing_values(df_pivot, selected_tasks):

    df_pivot = df_pivot.copy()
    imputed_rows = []

    for task in selected_tasks:
        missing_mask = df_pivot[task].isna()

        if missing_mask.any():
            missing_info = df_pivot.loc[missing_mask, ['Subject', 'Parcel']].copy()
            missing_info['Task'] = task
            imputed_rows.append(missing_info)

        parcel_medians = df_pivot.groupby('Parcel')[task].transform('median')
        global_median = df_pivot[task].median()

        df_pivot[task] = df_pivot[task].fillna(parcel_medians)
        df_pivot[task] = df_pivot[task].fillna(global_median)

    if len(imputed_rows) > 0:
        imputed_df = pd.concat(imputed_rows, ignore_index=True)

        imputed_summary = (
            imputed_df
            .groupby('Subject')
            .agg(
                n_imputati=('Parcel', 'count'),
                task_mancanti=('Task', lambda x: sorted(x.unique().tolist())),
                parcelle_mancanti=('Parcel', lambda x: sorted(x.unique().tolist()))
            )
            .reset_index()
        )
    else:
        imputed_df = pd.DataFrame(columns=['Subject', 'Parcel', 'Task'])
        imputed_summary = pd.DataFrame(
            columns=['Subject', 'n_imputati', 'task_mancanti', 'parcelle_mancanti']
        )

    return df_pivot, imputed_df, imputed_summary


def build_subject_task_matrix(df_pivot, subjects, parcels, selected_tasks):

    columns = []
    column_metadata = []

    for task in selected_tasks:
        for subj in subjects:
            subj_df = (
                df_pivot[df_pivot['Subject'] == subj]
                .set_index('Parcel')
                .loc[parcels]
            )

            values = subj_df[task].values.astype(float)

            columns.append(values)

            column_metadata.append({
                'column_index': len(column_metadata),
                'Subject': subj,
                'Task': task,
                'column_label': f"{task}_subj_{subj}"
            })

    X = np.column_stack(columns)
    column_metadata = pd.DataFrame(column_metadata)

    return X, column_metadata


def run_svd_on_subject_task_matrix(
    data_df,
    selected_tasks,
    output_dir,
    center_columns=True,
    impute_missing=True
):

    os.makedirs(output_dir, exist_ok=True)

    df_pivot, subjects, parcels = prepare_pivot_table(data_df, selected_tasks)

    if impute_missing:
        df_pivot, imputed_df, imputed_summary = impute_missing_values(
            df_pivot,
            selected_tasks
        )
    else:
        imputed_df = pd.DataFrame(columns=['Subject', 'Parcel', 'Task'])
        imputed_summary = pd.DataFrame(
            columns=['Subject', 'n_imputed', 'missing_tasks', 'missing_parcels']
        )

    X, column_metadata = build_subject_task_matrix(
        df_pivot=df_pivot,
        subjects=subjects,
        parcels=parcels,
        selected_tasks=selected_tasks
    )

    if np.isnan(X).any():
        raise ValueError(
            f"Matrix contains NaN values for tasks {selected_tasks}. "
        )

    if center_columns:
        X_centered = X - X.mean(axis=0, keepdims=True)
    else:
        X_centered = X.copy()

    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

    var_explained = (S ** 2) / np.sum(S ** 2)

    print(U.shape, S.shape, Vt.shape)

    # Save results
    np.savez(
        f"{output_dir}/svd_results.npz",
        U=U,
        S=S,
        Vt=Vt,
        var_explained=var_explained,
        X=X,
        X_centered=X_centered,
        parcels=np.array(parcels, dtype=object),
        subjects=np.array(subjects, dtype=object),
        tasks=np.array(selected_tasks, dtype=object)
    )

    column_metadata.to_csv(
        f"{output_dir}/column_metadata.csv",
        index=False
    )

    imputed_df.to_csv(
        f"{output_dir}/imputed_values_long.csv",
        index=False
    )

    imputed_summary.to_csv(
        f"{output_dir}/subjects_with_imputed_values.csv",
        index=False
    )

    summary = {
        'output_dir': output_dir,
        'selected_tasks': selected_tasks,
        'n_subjects': len(subjects),
        'n_parcels': len(parcels),
        'n_columns': X.shape[1],
        'matrix_shape': X.shape,
        'n_components': len(S),
        'variance_PC1': var_explained[0] if len(var_explained) > 0 else np.nan,
        'variance_PC2': var_explained[1] if len(var_explained) > 1 else np.nan,
        'n_subjects_with_imputed_values': len(imputed_summary)
    }

    return summary



def load_svd_result(config_dir, orient_components=False):

    npz_path = os.path.join(config_dir, "svd_results.npz")
    meta_path = os.path.join(config_dir, "column_metadata.csv")

    data = np.load(npz_path, allow_pickle=True)
    metadata = pd.read_csv(meta_path)

    U = data["U"].copy()
    S = data["S"].copy()
    Vt = data["Vt"].copy()
    var_explained = data["var_explained"].copy()
    tasks = data["tasks"].tolist()

    if orient_components:
        for k in range(min(U.shape[1], Vt.shape[0])):
            idx_max = np.argmax(np.abs(U[:, k]))
            if U[idx_max, k] < 0:
                U[:, k] *= -1
                Vt[k, :] *= -1

    return {
        "U": U,
        "S": S,
        "Vt": Vt,
        "var_explained": var_explained,
        "tasks": tasks,
        "metadata": metadata
    }


def load_svd_result_legacy(config_dir, orient_components=False):

    npz_path = os.path.join(config_dir, "svd_results.npz")
    meta_path = os.path.join(config_dir, "column_metadata.csv")

    data = np.load(npz_path, allow_pickle=True)
    metadata = pd.read_csv(meta_path)

    U = data["U"].copy()
    S = data["S"].copy()
    Vt = data["Vt"].copy()
    var_explained = data["var_explained"].copy()
    tasks = data["tasks"].tolist()

    if orient_components:
        for k in range(min(U.shape[1], Vt.shape[0])):
            idx_max = np.argmax(np.abs(U[:, k]))
            # if U[idx_max, k] < 0:
            #     U[:, k] *= -1
            #     Vt[k, :] *= -1

    return {
        "U": U,
        "S": S,
        "Vt": Vt,
        "var_explained": var_explained,
        "tasks": tasks,
        "metadata": metadata
    }


def compute_task_loading_summary(Vt, metadata, pc_index):

    row = Vt[pc_index, :]
    tmp = metadata.copy()
    tmp["loading"] = row

    grouped = tmp.groupby("Task")["loading"]

    summary = grouped.agg(["mean", "std"]).reset_index()
    summary["std"] = summary["std"].fillna(0.0)

    return summary, tmp
