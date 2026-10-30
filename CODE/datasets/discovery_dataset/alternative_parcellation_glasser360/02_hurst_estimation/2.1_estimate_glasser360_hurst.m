% Requires CODE/utils on MATLAB path (see README).
clear; clc;

%% Description
% This script estimates parcelwise Hurst exponents from fMRI data across multiple subjects
% and tasks.
% Parcellated data must be a .mat file of dimension length series x number
% of parcels contained in data_folder
%
% Requirements:
%   - nonfractal Toolbox (https://github.com/wonsang/nonfractal)
%   - wmtsa toolbox (https://www.atmos.washington.edu/~wmtsa/)
%   - MATLAB Statistics Toolbox

%%

addpath(genpath('wmtsa-matlab-0.2.6'));
addpath(genpath('nonfractal'));

%% parameters

% consts
data_folder = project_cluster_path('cluster_path_19'); % main folder containing data

% save params
verbose = true;

% bfn_mfin_ml params 
lb = [-0.5,0]; % lower bounds
ub = [1.5,10]; % upper bounds
filter = 'D4'; % wavelet transform filter
range = [1, 5];
wavelet = 'modwt';
padding = false;
decimate_b = false;

output_folder = sprintf('RESULTS/discovery_dataset/hurst/parcelwise/Glasser360/FIN_%s_%s_%d_%d_padding_%d_decimate_%d', filter, wavelet, range(1), range(2), padding, decimate_b); % path where data is saved

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%parallelization params
num_workers = 8;

%% 

% Find subjects automatically

subj_dirs = dir(fullfile(data_folder, 'sub-*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);

fprintf('Found %d subject folders in:\n%s\n\n', numel(subj_dirs), data_folder);

if isempty(subj_dirs)
    error('No subject folders found. Check data_folder.');
end

% start parallelization
if isempty(gcp('nocreate'))
           parpool('local', num_workers);
end

for j = 1:numel(subj_dirs)

    subject_id = subj_dirs(j).name;
    subject_folder = fullfile(data_folder, subject_id);

    fprintf('\n========================================\n');
    fprintf('Subject %d / %d: %s\n', j, numel(subj_dirs), subject_id);
    fprintf('========================================\n');

    % ---------------------------------------------------------------------
    % Output file
    % ---------------------------------------------------------------------
    H_filename = fullfile(output_folder, [subject_id '.mat']);

    if isfile(H_filename)
        fprintf('[SKIP] Hurst output already exists:\n%s\n', H_filename);
        continue;
    end

    % ---------------------------------------------------------------------
    % Find parcellated .mat file
    % ---------------------------------------------------------------------
    mat_files = dir(fullfile(subject_folder, '*w-noBPF-residuals*Glasser360*.mat'));

    if isempty(mat_files)
        % fallback: any .mat file in subject folder
        mat_files = dir(fullfile(subject_folder, '*.mat'));
    end

    if isempty(mat_files)
        warning('[SKIP] No parcellated .mat file found for %s in %s', ...
            subject_id, subject_folder);
        continue;
    end

    if numel(mat_files) > 1
        fprintf('Multiple .mat files found for %s. Using the first one:\n', subject_id);
        disp({mat_files.name}');
    end

    data_path = fullfile(mat_files(1).folder, mat_files(1).name);

    fprintf('Loading parcellated data:\n%s\n', data_path);

    % ---------------------------------------------------------------------
    % Load parcellated data
    % ---------------------------------------------------------------------
    S = load(data_path);

    if isfield(S, 'avg_time_series')
        data = S.avg_time_series;
    else
        % fallback: find first numeric 2D variable
        fn = fieldnames(S);
        data = [];

        for k = 1:numel(fn)
            candidate = S.(fn{k});
            if isnumeric(candidate) && ismatrix(candidate) && min(size(candidate)) > 1
                data = candidate;
                fprintf('Using variable "%s" as data.\n', fn{k});
                break;
            end
        end

        if isempty(data)
            warning('[SKIP] No valid numeric matrix found in %s', data_path);
            continue;
        end
    end

    [n_timepoints, n_parcels] = size(data);

    fprintf('Data size: %d timepoints x %d parcels\n', n_timepoints, n_parcels);

    Hursts = nan(1, n_parcels);

    % Store params 
    params.lb = lb;
    params.ub = ub;
    params.filter = filter;
    params.range = range;
    params.wavelet = wavelet;
    params.padding = padding;
    params.decimate_b = decimate_b;

    % ---------------------------------------------------------------------
    % Estimate Hurst exponent parcelwise
    % ---------------------------------------------------------------------
    parfor p = 1:n_parcels

        signal = data(:, p);

        if isempty(signal) || all(isnan(signal))
            Hursts(1, p) = NaN;
            continue;
        end

        % Handle partial NaNs if present
        if any(isnan(signal))
            signal = fillmissing(signal, 'linear', 'EndValues', 'nearest');
        end

        % Detrending and mean centering
        detrended_signal = detrend(signal, 'linear');
        detrended_signal = detrended_signal - mean(detrended_signal);

        if decimate_b
            detrended_signal = decimate(detrended_signal, 2);
        end

        T = length(detrended_signal);
        N_pad = 2^nextpow2(T);
        needed = N_pad - T;

        if needed > 0
            pad = flipud(detrended_signal(end-needed+1:end));
            signal_padded = [detrended_signal; pad];
        else
            signal_padded = detrended_signal;
        end

        if padding
            final_signal = signal_padded;
        else
            final_signal = detrended_signal;
        end

        try
            H_tmp = bfn_mfin_ml( ...
                final_signal, ...
                'range', range, ...
                'method', 'ML', ...
                'filter', filter, ...
                'lb', lb, ...
                'ub', ub, ...
                'verbose', 0, ...
                'wavelet', wavelet);

            Hursts(1, p) = H_tmp;

        catch ME
            Hursts(1, p) = NaN;
            fprintf('Error in %s parcel %d: %s\n', subject_id, p, ME.message);
        end

        if verbose
            fprintf('%s - Parcel %d / %d\n', subject_id, p, n_parcels);
        end
    end       
    
    save(H_filename, 'Hursts', 'subject_id', 'data_path', 'params');

end

delete(gcp('nocreate')); % close parallel pool if active



