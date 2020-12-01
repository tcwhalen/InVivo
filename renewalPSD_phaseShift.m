function [ psd_corr, phshift, psd_unc ] = renewalPSD_phaseShift(ts, op)
% Tim C Whalen, last edited July 2020
% Calculates renewal-corrected PSD and phase shift from spike train using
% Welch's method (average of multiple overlapping windows)
%
% Renewal-correction is inspired by the spike-shuffling method from Rivlin-
% Etzion et al 2006, "Local shuffling of spike trains boosts the accuracy 
% of spike train  spectral analysis", but uses an analytic (non-stochastic)
% improvement (see Gershner et al. "Neuronal Dynamics" textbook, ch.5
% https://neuronaldynamics.epfl.ch/online/Ch7.S5.html)
%
% Phase shift is the time-average of a coarse time-derivative of the phase 
% at each frequency in the input's fourier transform. It gives a measure of
% a frequency's "local stationarity", with zero indicating complete 
% stationary and pi indicating a wildly fluctuating signal. For details, 
% see Whalen et al 2020, "Delta oscillations are a robust biomarker of 
% dopamine depletion severity and motor dysfunction in awake mice".
%
% Due to edge effects, windows with fewer than 3 spikes are ignored.
%
% Outputs:
%   psd_corr: (:,1) double, renewal-corrected power spectral density
%       using Welch's method.
%   phshift: (:,1) double, mean phase shift (see Whalen 2020 methods)
%   psd_unc: (:,1) double, uncorrected power spectral density using
%       Welch's method)
%
% See arguments block for inputs. "op" are given as name-value pairs

% Bergman 
% Inputs:
%   ts: Nx1 num vector, spike times in sec
%   wind: int, window size for Welch's method (ideally power of 2, also
%       determines Rayleigh frequency)
%   step: int, step size for overlapping windows
%   FS: num, frequency to downsample spike train to (also determines
%       Nyquist frequency)

arguments
    %% Required arguments
    ts (:,1) double % spike times in seconds
    
    %% Optional tunable arguments
    op.wind double {mustBeInteger} = 2^12 % window size (in bins) for 
        % Welch's method (ideally power of 2, also determines Rayleigh 
        % frequency)
    op.step double {mustBeInteger} = 2^9 % step size (in bins) for 
        % overlapping windows. No overlap if step==wind
    op.FS double = 1000 % frequency (Hz) to downsample spike train to (also 
        % determines Nyquist frequency)
end

wind = op.wind; step = op.step; FS = op.FS;
freqs = 0:FS/wind:FS/2;

% convert to binary train
spkt = round(FS.*ts);
len_delt = spkt(end)-spkt(1)+1;
spkdelt = zeros(1,len_delt);
spkdelt( spkt - spkt(1)+1 ) = 1;

Ncalc = floor(len_delt/step - wind/step +1); % # fft's to do
psds = zeros(wind/2,Ncalc);
psds_unc = zeros(wind/2,Ncalc);
nspikes = zeros(1,Ncalc);
cphf = zeros(wind/2,Ncalc);

for s = 1:Ncalc
    segdelt = spkdelt(1+step*(s-1):step*(s-1)+wind);
    sps = find(segdelt);
    nspikes(s) = size(sps,2);
    if nspikes(s) <= 3 % edge case that gives anomalous results
        cphf(:,s) = NaN;
        psds(:,s) = NaN;
        psds_unc(s,:) = NaN;
        continue
    end
    
    % compute PSD of renewal process
    isi = diff(sps);
    edges = (0:wind)-.001;
    isin = histcounts(isi,edges);
    isidist = isin/sum(isin);
    phat = fft(isidist);
    phat(1)=0;
    chat = real((1+phat)./(1-phat));
    
    % compute PSD and normalize to renewal PSD
    ff = fft(segdelt-mean(segdelt));
    powf = abs(ff).^2; % averaging these ff's and scaling by 2/(FS*NFFT) gives exactly matlab's pwelch
    psd_corr = (powf/nspikes(s))./chat;
    psds_unc(:,s) = powf(1:wind/2)*2/(FS*wind);
    psds(:,s) = psd_corr(1:wind/2);
    
    % compute and align phases
    phf = angle(ff);
    phf = phf(1:size(psds_unc,1)); % remove repeated half;
    cphf(:,s) = mod(pi + (phf - 2*pi*(s-1)*step/1000*freqs(1:wind/2)),2*pi)-pi; % c=corrected, i.e. transformed so all start at zero
end

phdiff = diff(cphf,1,2); % phase differences over time
cphdiff = pi-abs(abs(phdiff)-pi); % make pi farthest, 2pi like zero
cphdiff = cphdiff(:,~isnan(cphdiff(1,:)));
phshift = mean(cphdiff,2);

nspikes_nonan = nspikes(~isnan(psds(2,:)));
nspikes_norm = nspikes_nonan/sum(nspikes_nonan);
psd_nonan = psds(:,~isnan(psds(2,:)));
psd_corr = sum(psd_nonan.*nspikes_norm,2); % rescale by #spikes in each window
psds_unc_nonan = psds_unc(:,~isnan(psds_unc(2,:)));
psd_unc = mean(psds_unc_nonan,2);
end