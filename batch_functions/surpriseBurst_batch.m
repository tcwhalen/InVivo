function data = surpriseBurst_batch(data)
% Tim Whalen, last edited Nov 2020
% Batch processing for surpriseBurst

% data must include: ts, T (in sec)
% To specify parameters, add the sub-struct "burst" to your data struct 
% containing any of the following you wish to define (all undefined will 
% use defaults):
%   fac = firing rate threshold (as multiple of baseline) to initially
%       consider burst (defualt = 1.5)
%   min_length_of_burst: minimum number of ISI's a burst can have
%       (default = 2)
%   local_length: time (in sec) in past over which to define current
%       baseline firing rate. If 0, calculates over entire train (default = 0)
%   surprise_cutoff: minimum Poisson surprise for a burst to be kept
%       (default = 10)

%% Initializations

% First get critical data; will throw error if not in struct
ts = data.ts;
T = data.T;

% Get parameters and defaults
names = {'fac', 'min_length_of_burst', 'local_length', 'surprise_cutoff'};
defaults = {1.5, 2, 0, 10}; % Wichmann
% defaults = {2, 4, 0, 2}; % Mastro et al. 2017
[data, fac, min_length_of_burst, local_length, surprise_cutoff] = extractDataInputs(data,'burst',names,defaults);

nfiles = length(ts);

if ~isfield(data,'burst')
    data.burst = struct();
end

data.burst.begin = cell(nfiles,1);
data.burst.num_spikes = cell(nfiles,1);
data.burst.surprise = cell(nfiles,1);
data.burst.rate = cell(nfiles,1);
data.burst.max_rate = cell(nfiles,1);
data.burst.baseline_rate = cell(nfiles,1);
data.burst.num_bursts = cell(nfiles,1);
data.burst.mean_spikes_per_burst = cell(nfiles,1);
data.burst.median_spikes_per_burst = cell(nfiles,1);
data.burst.total_spikes_in_bursts = cell(nfiles,1);
data.burst.mean_intra_burst_frequency = cell(nfiles,1);
data.burst.median_intra_burst_frequency = cell(nfiles,1);
data.burst.proportion_time_in_bursts = cell(nfiles,1);
data.burst.proportion_spikes_in_bursts = cell(nfiles,1);
data.burst.burst_duration = cell(nfiles,1);

%% Loop over files
for f = 1:nfiles
    disp(['Detecting bursts, file ' int2str(f) ' of ' int2str(nfiles)])
    nu = length(ts{f});
    
    data.burst.begin{f} = cell(nu,1);
    data.burst.num_spikes{f} = cell(nu,1);
    data.burst.surprise{f} = cell(nu,1);
    data.burst.rate{f} = cell(nu,1);
    data.burst.max_rate{f} = cell(nu,1);
    data.burst.baseline_rate{f} = cell(nu,1);
    data.burst.burst_duration{f} = cell(nu,1);
    
    data.burst.num_bursts{f} = zeros(nu,1);
    data.burst.bursts_per_second{f} = zeros(nu,1);
    data.burst.mean_spikes_per_burst{f} = zeros(nu,1);
    data.burst.median_spikes_per_burst{f} = zeros(nu,1);
    data.burst.total_spikes_in_bursts{f} = zeros(nu,1);
    data.burst.mean_intra_burst_frequency{f} = zeros(nu,1);
    data.burst.median_intra_burst_frequency{f} = zeros(nu,1);
    data.burst.proportion_time_in_bursts{f} = zeros(nu,1);
    data.burst.proportion_spikes_in_bursts{f} = zeros(nu,1);
    
    for u = 1:nu
        bs = surpriseBurst(ts{f}{u}, T(f), 'fac', fac, 'min_length_of_burst', min_length_of_burst, 'local_length', local_length, 'surprise_cutoff', surprise_cutoff);
        
        data.burst.begin{f}{u} = bs.begin;
        data.burst.num_spikes{f}{u} = bs.num_spikes;
        data.burst.surprise{f}{u} = bs.surprise;
        data.burst.rate{f}{u} = bs.rate;
        data.burst.max_rate{f}{u} = bs.max_rate;
        data.burst.baseline_rate{f}{u} = bs.baseline_rate;
        data.burst.burst_duration{f}{u} = bs.burst_duration;
        
        data.burst.num_bursts{f}(u) = bs.num_bursts;
        data.burst.bursts_per_second{f}(u) = bs.bursts_per_second;
        data.burst.mean_spikes_per_burst{f}(u) = bs.mean_spikes_per_burst;
        data.burst.median_spikes_per_burst{f}(u) = bs.median_spikes_per_burst;
        data.burst.total_spikes_in_bursts{f}(u) = bs.total_spikes_in_bursts;
        data.burst.mean_intra_burst_frequency{f}(u) = bs.mean_intra_burst_frequency;
        data.burst.median_intra_burst_frequency{f}(u) = bs.median_intra_burst_frequency;
        data.burst.proportion_time_in_bursts{f}(u) = bs.proportion_time_in_bursts;
        data.burst.proportion_spikes_in_bursts{f}(u) = bs.proportion_spikes_in_bursts;
    end
end