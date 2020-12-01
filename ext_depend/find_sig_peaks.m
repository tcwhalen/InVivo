function peak_inds = find_sig_peaks( data, thresh, max_n, srch_inds)
% function peak_inds = find_sig_peaks( data, thresh, local)
%
% Find local maxima in data that exceed threshold 
% Inputs:
%	data - vector containing data to search
%	thresh - threshold which data must exceed
%	max_n - identify local maxima over this number of points
%	srch_inds - indices of range w/in data over which to search
%
% Outputs:
%	peak_inds - indices into data of local maxima that exceed thresh
%
%	By RST 2005-07-31
%

max_inds = find_local_nmax( data, max_n);
sig_inds = find(data > thresh);
peak_inds = intersect( max_inds, sig_inds );

if exist('srch_inds','var')
	peak_inds = intersect( peak_inds, srch_inds);
end
