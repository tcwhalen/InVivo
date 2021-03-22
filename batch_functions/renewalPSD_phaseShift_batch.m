function [ data ] = renewalPSD_phaseShift_batch( data )
% Tim C Whalen, last edited July 2020
% Finds oscillations in spike trains and plots autocorrelation, renewal-
% corrected power spectral density (PSD) and mean phase shift plots.
%
% data must include: ts, T, animalcodees, rates, nfiles, files_trunc (for title, only if to_plot==1)
%
% To specify parameters, add the sub-struct "osc" to your data struct 
% containg any of:
%   acorr_bin: num, size of autocorrelation bin in sec (default = 0.02)
%   FS: num, frequency in Hz to downsample spike train to (also determines
%       Nyquist frequency and bins for step and wind) (default = 1000)
%   wind: int, window size in bins for Welch's method (ideally power of 2,
%       also determines Rayleigh frequency) (default = 2^12)
%   step: int, step size in bins for overlapping windows (default = 2^9)
%   srch_lo, hi: num, low and high values for range of frequencies in Hz 
%       from which to search for oscillation in Hz (default = 0.5, 4)
%   cntl_lo, hi: num, low and high values for range of frequencies in Hz
%       from which to derive confidence interval in Hz (high must be <=
%       FS/2) (default = 250, 500)
%   max_n: int, # points local exteme must be more extreme than to
%       count as peak/trough in PSD and phase shift plot (default = 7)
%   min_rate: num, skip neurons with firing rate lower than this (default = 5)
%   psd_threshp: num, p-value needed to flag significant oscillation, will
%   correct for multiple comparisons depending on size of srch_lo:srch_hi
%       (default = 0.01)
%   phase_threshp: num, p-value needed to flag significant oscillation, will
%   correct for multiple comparisons depending on number of PSD peaks found
%       (default = 0.05)
%   to_plot: bool, make acorr, PSD and phase shift plot for every unit
%       (default = 0)
%   examples: Nx2 int vector: file unit pairs to always plot, even if
%       to_plot==0 (default = [])
%
% The following fields will be added to osc for every unit using 
% the standard indexing structure {file}{unit} unless noted.
%   acorr: Nx1 num vector, spike autocorrelation
%   acorrtimes: Nx1 num vector, lags for plotting acorr
%   psd: Fx1 num vector, renewal-corrected power spectral density
%   psd_unc: Fx1 num vector, uncorrected (just Welch's) PSD
%   phase_shift: Fx1 num vector, mean phase shift
%   freqs: Fx1 num vector (not index, same for all units) frequency vector
%       for plotting psd, psd_unc, phase_shift
%   sigp_inds: Ux1 cell of ints, indices of freqs with significant 
%       oscillations determined by PSD and phase
%   sig_inds: Ux1 cell of ints, indices of freqs with significant 
%       oscillations determined by PSD only
%   has_osc: Ux1 bool vector, 1 if unit has osc determined by PSD and phase
%   has_osc_old: Ux1 bool vector, 1 if unit has osc determined only by PSD
%   frac_osc: Ax1 num vector (indexed by animal), fraction of oscillating
%       units for animal determiend by animalcodes
%   frac_osc_old: Ax1 num vector (indexed by animal), fraction of oscillating
%       units for animal determiend by animalcodes using only PSD

% get parameters or set to defaults
names = {'wind', 'step', 'acorr_bin', 'FS', 'srch_lo', 'srch_hi', 'cntl_lo', 'cntl_hi', 'min_rate', 'max_n', 'psd_threshp', 'phase_threshp', 'to_plot', 'examples'};
defaults = {2^12, 2^9, 0.02, 1000, 0.5, 4, 250, 500, 5, 7, .01, .05, 0, []};
[data, wind, step, acorr_bin, FS, srch_lo, srch_hi, cntl_lo, cntl_hi, min_rate, max_n, psd_threshp, phase_threshp, to_plot, examples] = extractDataInputs(data,'osc',names,defaults);

nfiles = data.nfiles;
sig_inds_all = cell(nfiles,1); 
sigp_inds_all = cell(nfiles,1); % both psd and phasediff passed
sig_freqs_all = cell(nfiles,1);
has_slow_osc = cell(nfiles,1);
has_slow_osc_old = cell(nfiles,1);
psd_corr_all = cell(nfiles,1);
psd_unc_all = cell(nfiles,1);
acorr_all = cell(nfiles,1);
phshift_all = cell(nfiles,1); % corrected phases

acorr_lag = 4000;

freqs = 0:FS/wind:50; % for most purposes only care about this range
freqs_long = 0:FS/wind:FS/2-FS/wind; % full range
data.osc.freqs = freqs;

% make set of examples complex for easy member checking
if ~isempty(examples)
    compex = complex(examples(:,1),examples(:,2));
else
    compex = [];
end

for f = 1:nfiles
    disp(['File ' int2str(f) ' of ' int2str(nfiles)])
    
    ts = data.ts{f};
    ts = ts(data.rates{f}>min_rate);
    rates = data.rates{f}(data.rates{f}>min_rate);
    
    nu = length(ts);
    sig_inds_all{f} = cell(nu,1);
    sigp_inds_all{f} = cell(nu,1);
    sig_freqs_all{f} = cell(nu,1);
    psd_corr_all{f} = zeros(nu,length(freqs));
    psd_unc_all{f} = zeros(nu,length(freqs));
    phshift_all{f} = zeros(nu,length(freqs));
    has_slow_osc{f} = zeros(nu,1);
    has_slow_osc_old{f} = zeros(nu,1);
    acorr_all{f} = zeros(nu,acorr_lag/(acorr_bin*1000)+1);
    
    for u = 1:length(ts)
        % Compute autocorr
        [acorr,acorrlag] = spikeAutocorr(ts{u},acorr_lag, acorr_bin);
        acorr_all{f}(u,:) = acorr;
        
        % corr is renewal-corrected, unc is uncorrected (basic Welch's)
        [psd_corr, phshift, psd_unc] = renewalPSD_phaseShift(ts{u}, 'wind', wind, 'step', step, 'FS', FS);
        
        srch_inds = [find(freqs_long>=srch_lo,1), find(freqs_long<=srch_hi,1,'last')];
        cntl_inds = [find(freqs_long>cntl_lo,1), find(freqs_long<=cntl_hi,1,'last')];
        [sigp_inds,sig_inds] = findSigOsc(psd_corr,phshift,srch_inds,cntl_inds,max_n,psd_threshp,phase_threshp);

        if to_plot || ismember(complex(f,u),compex)
            figure
            subplot(3,1,1)
            plot(acorrlag,acorr(1:length(acorrlag)),'k','LineWidth',2)
            xlabel('Time (ms)')
            ylabel('Acorr.')
            makeNice(gca)
            title('Autocorrelation')
            
            subplot(3,1,2)
            hold on
            plot(freqs_long,psd_corr(1:length(freqs_long)),'k','LineWidth',2)
            plot(freqs(sigp_inds), psd_corr(sigp_inds),'ro','LineWidth',2);
            xlim([.1 16])
            xlabel('Frequency (Hz)')
            ylabel('Power')
            title('Renewal-Corrected PSD')
            makeNice(gca)
            
            subplot(3,1,3)
            hold on
            plot(freqs_long(2:end),phshift(2:length(freqs_long),:),'k','LineWidth',2)
            plot(freqs(sigp_inds), phshift(sigp_inds),'ro','LineWidth',2);
            xlim([.1 16])
            xlabel('Frequency (Hz)')
            ylabel('\Deltaphase')
            title('Phase Shift')
            makeNice(gca)
        end
        sig_freqs_all{f}{u} = freqs(sig_inds);
        sig_inds_all{f}{u} = sig_inds;
        sigp_inds_all{f}{u} = sigp_inds;
        psd_corr_all{f}(u,:) = psd_corr(1:length(freqs));
        psd_unc_all{f}(u,:) = psd_unc(1:length(freqs));
        phshift_all{f}(u,:) = phshift(1:length(freqs));
        if ~isempty(sigp_inds)
            has_slow_osc{f}(u) = 1;
        end
        if ~isempty(sig_inds)
            has_slow_osc_old{f}(u) = 1;
        end
    end
end

nanimals = max(data.animalcodes);
frac_osc = zeros(nanimals,1);
frac_osc_old = zeros(nanimals,1);
for i = 1:nanimals
    osc = cellfun(@(x) ~isempty(x), getit(sigp_inds_all(data.animalcodes==i)));
    osc_old = cellfun(@(x) ~isempty(x), getit(sig_inds_all(data.animalcodes==i)));
    frac_osc(i) = mean(osc);
    frac_osc_old(i) = mean(osc_old);
end

data.osc.sigp_inds = sigp_inds_all;
data.osc.sig_inds = sig_inds_all;
data.osc.has_osc = has_slow_osc;
data.osc.has_osc_old = has_slow_osc_old;
data.osc.frac_osc = frac_osc;
data.osc.frac_osc_old = frac_osc_old;
data.osc.acorr = acorr_all;
data.osc.acorrtimes = acorrlag;
data.osc.psd = psd_corr_all;
data.osc.psd_unc = psd_unc_all;
data.osc.phase_shift = phshift_all;
data.osc.freqs = freqs;