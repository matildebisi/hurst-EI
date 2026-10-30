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
addpath('nonfractal');

%% parameters

% consts
subj_ids = 1:36; % subject ids in file names
tasks = {'00', '01', '03', '05'};
data_folder = 'DATA/parcellated_residuals'; % main folder containing data

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

output_folder = sprintf('RESULTS/Hurst_results/residuals/parcelwise/FIN_%s_%s_%d_%d_padding_%d_decimate_%d', filter, wavelet, range(1), range(2), padding, decimate_b); % path where data is saved

%parallelization params
num_workers = 6;

%% 

% start parallelization
if isempty(gcp('nocreate'))
           parpool('local', num_workers);
end

for j = 1:numel(subj_ids) % for each subject
    id = subj_ids(j);   

    for i = 1:length(tasks) % for each task
        data_path = sprintf('%s/sub%02d/task_%s.mat', data_folder, id, tasks{i}); % parcellated data path
        
        if isempty(data_path) || ~isfile(data_path)
            warning('File not found: %s', data_path);
            continue;
        end
        
        % load parcellated data file
        S = load(data_path);
        fn = fieldnames(S);
        data = S.(fn{1});
        [l, n_parcels] = size(data); 

        Hursts = zeros(1, n_parcels);
        likelihoods = zeros(1, n_parcels);
        
        % Estimate Hurst exponent using maximum likelihood method
        parfor p = 1:n_parcels
            signal = data(:, p); 
            
            % detrending and mean centering
            detrended_signal = detrend(signal, 'linear'); % removes linear trend
            detrended_signal = detrended_signal - mean(detrended_signal); % mean-centering

            if decimate_b
                detrended_signal = decimate(detrended_signal, 2)
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

            if isempty(detrended_signal)
                Hursts(1, p) = NaN;
                likelihoods(1, p) = NaN;
                continue;
            end

            if padding
                final_signal = signal_padded;
            else
                final_signal = detrended_signal;
            end

            [Hursts(1, p), ~, ~, ~] = bfn_mfin_ml(final_signal, 'range', range, 'method', 'ML','filter', ...
                filter, 'lb', lb, 'ub', ub, 'verbose', 0, 'wavelet', wavelet); 

            if verbose
                fprintf('Sub %d - Parcel %d / %d\n', id, p, n_parcels);
            end
        end
        % save data in output_folder
        H_filename = sprintf('%s/sub%02d_task%s.mat', output_folder, id, tasks{i});

        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end

        save(H_filename, 'Hursts');
    end
    
end

delete(gcp('nocreate')); % close parallel pool if active



