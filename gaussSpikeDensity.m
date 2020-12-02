function [ times, sdf ] = gaussSpikeDensity( ts, std, lo, hi, op)
% Tim Whalen 7/3/17, last edited Nov 2020
% Convolves a spike train (impulses) with a Gaussian to generate the spike
% density function (i.e. smoothed instantaneous firing rate)
%
% Output
% times: x-axis to plot result against (sec)
% sdf: spike density function
%
% See arguments block for inputs. "op" are given as name-value pairs

arguments
    % Required arguments
    ts (:,1) double % spike times in seconds
    std double % standard deviation of Gaussian smoothing kernel
    lo double % time in sec to begin smoothing. Only information after this
        % time is used, and the first time returned will be 3*std + lo to
        % avoid edge effects while maintaining no leakage from before lo
    hi double % time in sec to end smoothing. Same details apply as for lo; 
        % no information after hi is used and 3*std is cut
        
    % Optional tunable arguments
    op.FS double = 1000 % sampling frequency of returned time series in Hz 
        % (i.e. 1/bin length) Should be less than the aquisition FS of ts.
    op.cut double = 3 % number of std's to remove from either end of signal
        % after smoothing to avoid edge effects. Default of 3 gets you
        % 99.7% of the way there; any lower gives a longer result if your
        % input ts is very short, but interpret edges at your own risk.
end

FS = op.FS; cut = op.cut;

dt = 1/FS; % downsample to 1 kHz
times_full = lo:dt:hi;
nbins = length(times_full);
binm = 1/dt; % bins/sec

delt = zeros(1,nbins);
delt( round(binm*ts(ts<hi))+1) = 1;
gauss = normpdf(-std*cut:dt:std*cut,0,std);

sdf_full = conv(delt,gauss,'same'); % same keeps only length of A
sdf = sdf_full(binm*std*cut+1:end-binm*std*cut);
times = (lo+std*cut):dt:(hi-std*cut);
end