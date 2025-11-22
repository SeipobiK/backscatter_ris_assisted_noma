%% pre_generate_channels.m
clear all; clc;
addpath(genpath('/home/morolong/Documents/MATLAB/backscatter_ris_assisted_noma'));
rng(2023, 'twister');


% Initialize parameters
para = para_init();
[BS_array, RIS_array] = generate_arrays(para);

MC_MAX = para.MC_MAX;
rng_seeds = randi(1e6, MC_MAX, 1); % Seeds for reproducibility

% Preallocate channel containers
H_all = cell(MC_MAX,1);
g_all = cell(MC_MAX,1);
f_all = cell(MC_MAX,1);

% Generate channels
for mc = 1:MC_MAX
    
    rng(rng_seeds(mc), 'twister');
    [H_all{mc}, g_all{mc}, f_all{mc}] = generate_channel(para, BS_array, RIS_array);
end

% Save channels to MAT file
% Create a folder for storing pre-generated channels
channel_dir = 'pre_generated_channels';
if ~exist(channel_dir, 'dir')
    mkdir(channel_dir);
end

% Construct a clean, descriptive filename
timestamp = datestr(now,'yyyymmdd_HHMMSS');
filename = fullfile(channel_dir, ...
    sprintf('channels_M%d_N%d_%s.mat', para.MC_MAX, para.N, timestamp));

% Save all channel data
save(filename, 'H_all', 'g_all', 'f_all', 'rng_seeds', 'para', 'BS_array', 'RIS_array');

fprintf('Pre-generated channels saved to: %s\n', filename);
