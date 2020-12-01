function [ ecog_reg ] = ecogSpikeRegress(ts_cell, ecog, T, ecog_start, ecog_FS, op)
% Computes autoregressive model for ECoG (null) and adds in parameter for
% several delays of a Gaussian smoothed spike train (spike density 
% function) to find maximum likelihood time lag (or lead) between ECoG and 
% neuron and if it is significantly better than a null autoregressive
% model. Computes on the difference of ECoG samples (i.e. ARI model with d=1)
%
% Takes a cell array of spike times from simutlaneously recorded units, so
% it can process an entire file at once only building one null AR model of
% the ECoG signal.
%
% Note that a new model is computed for each lag separately - this is not a
% distirbuted lag model, as multicolinearity with the SDF would heavily
% bias regression coefficients. Building a separate model for each lag is
% reasonable under the assumption that if a lag exists by which the neuron
% firing influences the ECoG or vice versa, then there is only one such lag
% by which this influence occurs.
%
% See arguments block for input. "op" are given as  name-value pairs
%
% Returns a struct witht he following fields (where nu is the number of
% units in ts_cell and nlag is nlag_past+nlag_fut+1)
% sig_lags: double (nu,1) time of best lag where negative is spike
%   predicting ECoG, positiv is vice-versa. NaN if no significance
% bs: double (nu,nlag): regression coefficients
% ssrs: double (nu,nlag): sum of squared residuals
% rmse_scaled: double (nu,nlag): root mean squared error scaled to mean
%   ecog emplitude
% pvals: double (nu,1): p-value for best lag, corrected for multiple
%   comparisons (so may be >1)
% pttype: int (nu,1): 1 if unit's best lag is at a peak, -1 if trough, 0 if
%   sig_lag = 0 for this unit
% mean_amp: double: mean absolute amplitude of zero-centered ecog

arguments
    %% required
    ts_cell (:,1) cell % each element is a (:,1) vector of spike times in 
        % sec from different simultaneously recorded units. If using a
        % single unit, still wrap the vector in a {1,1} cell
    ecog (:,1) double % ecog time series. NaN's are acceptable for missing/
        % noisy data
    T double % length of recording in sec
    ecog_start double % time in sec of first bin of ecog
    ecog_FS double % ECoG sampling frequency
    
    %% tunable arguments
    op.bin double = 0.1 % bin width in sec
    op.step double = 0.1 % time difference between bins (set equal to bin for 
        % no overlap)
    op.ar_lag double {mustBeInteger} = 1 % integer number of terms in AR model
    op.nlag_past double {mustBeInteger} = 100 % number of lags (in units of
        % bins) to make alternate models from using past spike information.
        % Must be > ar_lag
    op.nlag_fut double {mustBeInteger} = 100 % number of lags (in units of
        % bins) to make alternate models from using future spike information
    op.alpha double = 0.05 % significance threshold for comparison to AR
        % model, which will beautomatically corrected for multiple 
        % comparisons (i.e. *1/(nlag_past + nlag_fut + 1)
    op.do_sdf logical = 1 % use Gaussian-smoothed spike trains instead of point
        % processes. 
    op.sdf_std double = 0.1 % standard deviation of Gaussian kernel for sdf
        % computation in seconds
    op.backwards logical = 0 % 1 to train AR model in backwards time
end

bin = op.bin; step = op.step; ar_lag = op.ar_lag; nlag_past = op.nlag_past;
nlag_fut = op.nlag_fut; alpha = op.alpha; do_sdf = op.do_sdf; sdf_std = op.sdf_std;
backwards = op.backwards;

if nlag_past< ar_lag
    error('need nlag_pos > ar_lag')
end
nlag_total = nlag_past+nlag_fut+1; % +1 for lag zero

%% Downsample ecog to bin widths
if do_sdf
    sdf_FS = 1000;
    binstarts = (3*sdf_std:step:T-3*sdf_std-bin)';
else
    binstarts = (ecog_start:step:T-bin-.1)'; % .1 fixes lost signal at end
end

ecogrs = zeros(size(binstarts)); % rs = resampled
binstartsind = floor((binstarts-ecog_start)*ecog_FS+1); % indices for ecog
binind = floor(bin*ecog_FS+1);
for i = 1:length(binstarts)
    ecogrs(i) = mean(ecog(binstartsind(i):binstartsind(i)+binind));
end
ecogdiff = [0; diff(ecogrs)];
if backwards
    ecogdiff = flipud(ecogdiff);
end

%% Construct autoregressive model of ECoG
npoints = length(ecogdiff)-nlag_total+1; % # data points in regression - remove enough for maxlag once spk data is added
nparams = ar_lag+1;
y = ecogdiff(nlag_past+1:end-nlag_fut);
ar_data = zeros(npoints,ar_lag+1);
ar_data(:,1) = ones(npoints,1); % intercept
for i = 1:ar_lag % unintuitive but faster way to fill in regressor matrix
    ar_data(:,i+1) = ecogdiff(i+nlag_past-ar_lag:end-ar_lag+i-1-nlag_fut);
end

% regress automatically ignores NaN rows, but we'll remove them anyway
badrowsa = any(isnan(ar_data),2);
badrowsy = isnan(y);
badrows = (badrowsa+badrowsy)>=1; % union
ar_data(badrows,:)=[];
y(badrows,:)=[];

[~, ~, r_null] = regress(y,ar_data);
ssr_null  = r_null'*r_null; % sum of squared residuals

%% Construct models each with different single spike train lag
nu = size(ts_cell,1);
sig_lags = zeros(length(ts_cell),1)+NaN;
bs_all = zeros(length(ts_cell),nlag_total);
ssrs_all = zeros(length(ts_cell),nlag_total);
rmse_scaled_all = zeros(length(ts_cell),nlag_total);
pvals = zeros(length(ts_cell),1);
pttype = zeros(length(ts_cell),1);

mean_amp = mean(abs(ecog(~isnan(ecog))));
for u = 1:nu
    ts = ts_cell{u};
    ssrs = zeros(nlag_total,1);
    bs = zeros(nlag_total,ar_lag+2); % coefficients
    
    if do_sdf
        [sdftimes, sdf] = gaussSpikeDensity( ts, sdf_std, 0, T);
        if sdftimes(1)~= sdf_std*3
            error('sdftimes(1)~= sdf_std*3')
        end
        binstartsind_sdf = floor((binstarts-sdf_std*3)*sdf_FS+1.1); % 1000 hard-coded into gaussSpikeDensity. Extra 0.1 (1.1 vs just 1) fixes rounding errors
    end
    
    tsrs = zeros(size(binstarts));
    for i = 1:length(binstarts)
        if do_sdf
            tsrs(i) = mean(sdf(binstartsind_sdf(i):binstartsind_sdf(i)+bin*sdf_FS-1));
        else
            tsrs(i) = sum(ts>=binstarts(i) & ts<binstarts(i)+bin);
        end
    end
    if backwards
        tsrs = flipud(tsrs);
    end
    
    lags = -nlag_past:nlag_fut;
    for i = 1:nlag_total
        lag = lags(i); % as lag increases, going FORWARDS in time
        x = tsrs(nlag_past+1+lag:end+lag-nlag_fut);
        x(badrows)=[];
        [bs(i,:), ~, r] = regress(y,[ar_data x-mean(x)]);
        ssrs(i) = r'*r; % sum of squared residuals
    end
    
    [ssr_best, ssr_besti] = min(ssrs);
    ssr_best_time = (ssr_besti-nlag_past-1)*bin; % convert from index to time
    
    Fstat = (ssr_null-ssr_best)/(ssr_best/(npoints-nparams)); % simplified because null has exactly 1 fewer parameter
    Fpval = 1-fcdf(Fstat,1,npoints-nparams);
    pttype(u) = sign(bs(ssr_besti,end))*(Fpval<alpha/nlag_total);
    
    bs_all(u,:) = bs(:,end);
    ssrs_all(u,:) = ssrs;
    rmse_scaled_all(u,:) = sqrt(ssrs/npoints)/mean_amp;
    pvals(u) = Fpval*nlag_total;
    if Fpval<alpha/nlag_total
        sig_lags(u) = ssr_best_time;
    end
end

ecog_reg = struct();
ecog_reg.sig_lags = sig_lags;
ecog_reg.bs = bs_all;
ecog_reg.ssrs = ssrs_all;
ecog_reg.rmse_scaled = rmse_scaled_all;
ecog_reg.pvals = pvals;
ecog_reg.pttype = pttype;
ecog_reg.mean_amp = mean_amp;
end