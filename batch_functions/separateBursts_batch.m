function [ data ] = separateBursts_batch( data )
% Tim Whalen, last edited Nov 2020
% Makes new spike time vectors removing or only including bursts
% data must include: ts, T, nfiles, burst.begin, burst.num_spikes
% burst struct can be created by running surpriseBurst
% No parameters necessary

% Output: The following fields will be added to data.burst
%   ts_burst: spike times only for spikes in bursts
%   ts_noburst: spike times with all bursts removed
%   bin_burst: binary vector with 1's during bursts and 0's elsewhere,
%     sampled at 1 kHz (i.e., plot against 0:1/1000:T

nfiles = data.nfiles;
tsall = data.ts;
Tall = data.T;
startsall = data.burst.begin;
nspikesall = data.burst.num_spikes;

dt = 1/1000;

ts_burst = cell(nfiles,1);
ts_noburst = cell(nfiles,1);
bin_burst = cell(nfiles,1);

for f = 1:nfiles
    ts_burst{f} = cell(size(tsall{f}));
    ts_noburst{f} = cell(size(tsall{f}));
    bin_burst{f} = cell(size(tsall{f}));
    for u = 1:length(tsall{f})
        ts = tsall{f}{u};
        T = Tall(f);
        starts = startsall{f}{u};
        nspikes = nspikesall{f}{u};
        
        bi = zeros(1,sum(nspikes)); % indices of burst spikes
        bin_burst{f}{u} = zeros(length(0:dt:T),1);
        count = 1;
        for i = 1:length(starts)
            bi(count:count+nspikes(i)-1) = starts(i):starts(i)+nspikes(i)-1;
            bin_burst{f}{u}(round(ts(starts(i))/dt)+1:round(ts(starts(i)+nspikes(i)-1)/dt)) = 1;
            count = count+nspikes(i);
        end
        ts_burst{f}{u} = ts(bi);
        ts_noburst{f}{u} = ts(setdiff(1:length(ts),bi));
    end
end

data.burst.ts_burst = ts_burst;
data.burst.ts_noburst = ts_noburst;
data.burst.bin_burst = bin_burst;

end