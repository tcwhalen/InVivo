function [sigp_inds,sig_inds] = findSigOsc(psd,phshift,srch_inds,cntl_inds,max_n,psd_threshp,phase_threshp)
% Tim C Whalen, last edited July 2020
% Given PSD and phase shift, finds significant oscillations
% Inputs:
%   psd: Nx1 vector, power spectrum (N = # frequencies)
%   phshift: Nx1 vector, phase shift plot
%   srch_inds: 1x2 int, low and high indices of psd/phshift to search for 
%       significant peak/trough
%   cntl_inds: 1x2 int, low and high indices to use for control confidence 
%       interval
%   max_n: int, # points local exteme must be more extreme than to
%       count as peak/trough
%   psd, phase_threshp: num, p-value threshold for significance
% Outputs:
%   sigp_inds: Nx1 int vector, indices of significant phshift trough and 
%       psd peak
%   sig_inds: Nx1 int vector, indices of significant psd peak, ignoring 
%       phshift

p = psd_threshp/(srch_inds(2)-srch_inds(1)+1);
thresh = norminv(1-p) * std( psd(cntl_inds(1):cntl_inds(2)) ) + 1;
sig_inds = find_sig_peaks(psd(1:srch_inds(end)+ceil(max_n/2)), thresh, max_n, srch_inds(1):srch_inds(2));

phasep = phase_threshp/size(sig_inds,1);
phase_thresh = -norminv(1-phasep) * std(phshift(cntl_inds(1):cntl_inds(2)) ) + mean(phshift(cntl_inds(1):cntl_inds(2)));
sigp_inds = sig_inds(phshift(sig_inds)<phase_thresh);
