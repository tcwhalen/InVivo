function [data] = spikeSpectrogram(ts, op)
% Tim Whalen, last edited Nov 2020
% Returns and plots renewal-corrected spectrogram of spike train
% Will optionally return an additional 1D time series "series" (e.g. 
% simultaneously recorded LFP or mvmt) binned at same times as spectrogram
%
% Returns a struct with the following fields:
% spec (:,:): the spectrogram of ts. spec(f,t) gives the power for the f-th
%   frequency at the t-th time
% spec_smooth (:,:): same as spec, smoothed with a 3x3 Gaussian filter
% times (:,1): the time axis (center of bin) for spec in seconds
% freqs (:,1): the frequency axis for spec in Hz
% pair_means (:,1): only if pair_sig included, pair_sig binned at "times"
%
% See arguments block for input. "op" are given as name-value pairs
% Exampelc all: spikeSpectrogram(mySpikes, 'FS', 5000, 'FR_norm', 1)
% Plot output with imagesc(times, freqs, spec_smooth)
% Recommended to call set(gca,'ydir','normal') to flip y-axis

arguments
    %% required inputs
    ts (:,1) double % spike times, in sec
    
    %% tunable arguments
    op.wind double {mustBeInteger} = 2^13; % number of bins in FFT window
    op.step double {mustBeInteger} = 2^11; % step size for overlapping
        % windows in bins, no overlap if wind==step
    op.FS double = 1000; % Hz to downsample ts to
    op.HF_norm double = 200 % normalize each bin to the mean of the high 
        % frequency band from this value (Hz) to the Nyquist frequency. If
        % zero, this correction isn't used. Must be < FS/2, but too close
        % to FS/2 will use few samples and give unstable results
    op.FR_norm logical = 0 % normalize each bin by its local firing rate (not
        % to be confused with renewal-correction, which is always used)
        % Using alongside FR_norm is not recommended
    
    %% optional: return a paired signal binned alongside spectrogram
    op.pair_sig (:,1) double % optional paired signal to downsample to
        % same bins as spectrogram. Must be time series, not timestamps
        % like ts. Ignored if empty
    op.pair_start double % time of first pair_sig bin
    op.pair_FS double % sampling frequency of pair_sig
    op.pair_norm double = 0 % 1 to normalize pair_sig by z-score
end

do_pair = 0;
if isfield(op,'pair_sig')
    if ~isfield(op,'pair_FS') || ~isfield(op,'pair_start')
        error('if specifying pair_sig, must also specify pair_start and pair_FS')
    end
    do_pair = 1;
    pair_sig = op.pair_sig; pair_start = op.pair_start;
    pair_FS = op.pair_FS; pair_norm = op.pair_norm;
    pair_t = (pair_start:1/pair_FS:length(pair_sig))';
end
wind = op.wind; step = op.step; FS = op.FS; HF_norm = op.HF_norm; FR_norm = op.FR_norm;

spkt = round(FS.*ts);
len_delt = spkt(end)-spkt(1)+1;
spkdelt = zeros(1,len_delt);
spkdelt( spkt - spkt(1)+1 ) = 1;

Ncalc = floor(len_delt/step - wind/step +1); % # fft's to do
spec = zeros(wind/2,Ncalc);
full_freqs = 0:FS/wind:FS/2;
if HF_norm
    min_norm = find(full_freqs>=FR_norm,1);
end
windfun = hamming(wind)';
if do_pair
    pair_means = zeros(Ncalc,1);
end

for s = 1:Ncalc
    if do_pair
        pair_seg = abs(pair_sig(pair_t>=(1+step*(s-1))/FS & pair_t<(step*(s-1)+wind)/FS));
        pair_seg(pair_seg==1)=0;
        pair_means(s) = mean(pair_seg);
    end
    
    deltwin = spkdelt(1+step*(s-1):step*(s-1)+wind);
    
    % compute PSD of renewal process
    isi = diff(find(deltwin));
    edges = (0:wind)-.001;
    isin = histcounts(isi,edges);
    isidist = isin/sum(isin);
    phat = fft(isidist);
    chat = real((1+phat)./(1-phat));
    
    % compute PSD and normalize to renewal PSD
    ff = abs(fft(windfun.*(deltwin-mean(deltwin)))).^2/(FS/2);
    psd_corr = ff./chat;
    spec(:,s) = psd_corr(1:wind/2);
    if HF_norm
        spec(:,s) = spec(:,s)./mean(spec(min_norm:end,s));
    end
    if FR_norm
        spec(:,s) = spec(:,s)/sum(deltwin);
    end
end
if do_pair && pair_norm
    pair_means = zscore(pair_means);
end

smooth2 = fspecial('gaussian', [3 3], .6); % 3x3 Gaussian kernel
smooth2 = smooth2/sum(sum(smooth2));
spec_smooth = conv2(spec,smooth2,'same');
        
cents = (wind/2:step:wind/2+step*(Ncalc-1))/1000; % times, at center of bins

data.spec = spec;
data.spec_smooth = spec_smooth;
data.times = cents;
data.freqs = full_freqs;
if do_pair
    data.pair_means = pair_means;
end