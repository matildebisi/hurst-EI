clear; clc; close all;

%% Description
% This script estimates parcel-wise Hurst exponents from fMRI data across multiple subjects
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
tasks = {'01', '02', '03', '04', '05', '06'};
label = '';
data_folder = 'DATA/pupillometry/pupil_residuals/pupil_size/decimated_2s'; % main folder containing data

% bfn_mfin_ml params 
lb = [-0.5,0]; % lower bounds
ub = [1.5,0]; % upper bounds
maxit = 100;
range = [1, 5];
filter = 'D4'; % wavelet transform filter
omegamode = 'cov';
wavelet = 'modwt';

% parallelization params
parallelize = false;
num_workers = 12;

% save params
verbose = true;
output_folder = 'RESULTS/replication_dataset/hurst/eye_tracking/pupil_size/decimated_2s'; % path where data is saved
H_filename = sprintf('%s/FIN_hurst_exponents_%s_%s_%d_%d.mat', output_folder, filter, wavelet, range(1), range(2));

%% 
% figure;
% hold on;

% start parallelization
if parallelize
    if isempty(gcp('nocreate'))
                parpool('local', num_workers);
    end
end

hursts = ones(length(subj_ids), length(tasks))*(-1);

for j = 1:numel(subj_ids) % for each subject
    id = subj_ids(j);   

    for i = 1:length(tasks) % for each task
        data_path = sprintf('%s/sub%02d_%s%s.mat', data_folder, id, tasks{i}, label); % data path
        
        if isempty(data_path) || ~isfile(data_path)
            warning('File not found: %s', data_path);
            continue;
        end
        
        % load parcellated data file
        S = load(data_path);
        fn = fieldnames(S);
        data = S.(fn{1});
%         N  = length(data);
%         N2 = 2^floor(log2(N));
%         data = data(end-N2+1:end);

      %  plot(data)
        
        input = double(data');
        [hursts(j, i), ~, ~, ~] = bfn_mfin_ml(input, 'method', 'ML','filter', ...
            filter, 'lb', lb, 'ub', ub, 'verbose', 0, 'maxit', maxit, 'wavelet', wavelet, ...
            'omegamode', omegamode, 'd_only', true, 'abstol', 10e-14, 'range', range); 

        if verbose
            fprintf('Sub %d - Task %d \n', id, i);
        end
    end
end
   
% save data in output_folder

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
save(H_filename, 'hursts');

delete(gcp('nocreate')); % close parallel pool if active



