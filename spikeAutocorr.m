function [y,t] = spikeAutocorr(A,lag,binsize,datalen)
% 06-04-14 by RST
% Edited Tim Whalen 1/9/17 to include variable binsize (lag still in sec)
%
%	Computes the autocorrelation
%	of input spike times A (using lag samples).
%	The output t is the lag time.
%
% Sample call:
% 		[B t] = autocorr(A,1000,0.001,20);
%
% Inputs:
%	A - spike times in seconds
%	lag - number of msec bins (default = 500)
%   binsize - width of bins in sec (default = 0.001)
%	datalen - length (in sec) of A

if ~exist('lag','var')
	lag = 500;
end
if ~exist('datalen','var')
	datalen = max(A);
end
if ~exist('binsize','var')
	binsize = .001; % 1msec bins
end
FS = 1/binsize;
A = A(:);

% Make delta functions
deltlen = round(FS*datalen);
delt = zeros(1,deltlen);
delt( round(FS*A)+1) = 1;

% calculation of autocorrelation
binlag = lag/(binsize*1000); % converts msec lag to #bins
[y,t] = xcorr(delt,delt,binlag);
t = t*binsize*1000; % convert back to msec

y = FS * y ./ length(A);			% calibrate to spike/sec

y = y(binlag+1:2*binlag+1);	% removal of negative half of autocorrelation
y(1) = 0;
t = t(binlag+1:2*binlag+1);