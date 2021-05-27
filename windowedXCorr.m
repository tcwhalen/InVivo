function [ xco_avg, xcos ] = windowedXCorr(sig1, sig2, s_len, b_len, op)
% Tim Whalen, last edited May 2021
% Computes a modified cross-correlation as an average of xcorrs computed
% over a moving window
%
% In brief, the algorithm is:
% 1. Window sig1 and sig2 to length b_len starting at time zero
% 2. Mean-subtract the windowed signals
% 3. Zero sides of windowed sig2 so only the middle s_len values remain
% 4. Compute xcorr with maximum lag = (b_len-s_len)/2
% 5. Move the original window over by length step_size for both signals
% 6. Repeat 2-5 until end of signals are reached
% 7. Average all xcorr's to acheive final result
%
% The moving window corrects for nonstationarities at timescales larger
% than b_len. The zeroing in each window corrects for the artificial
% decrease in variance of a cross-correlation at increasing lags (maxlag is
% chosen such that an equal number of zeroed valeus are used at every lag)
%
% See arguments block for input. "op" are given as name-value pairs
%
% Returns:
% xco_avg: (1,:) double, cross-correlation averaged across windows. Length
%   is b_len-s_len+1
% xcos: (:,:) double, cross-correlation for every window. First dim is
%   window index, second is lag

arguments
    %% required
    sig1 (:,1) double % signal #1 for cross-correlation
    sig2 (:,1) double % signal #2 for cross-correlation. Should be same
        % size and sampling frequency as sig1
    s_len double {mustBeInteger} % length of signal to keep unzeroed in
        % each window
    b_len double {mustBeInteger} % window length. maxlag is set to
        % (b_len-s_len)/2
        
    %% optional
    op.step_size double {mustBeInteger} = b_len % length to move window by 
        % each step. Default is no overlap
    op.to_plot logical = 0 % 1 to plot resulting xcorr
    op.FS double = 1000 % samplign frequency for x-axis of plot. Default
        % assumes 1 kHz signal
end
step_size = op.step_size; to_plot = op.to_plot; FS = op.FS;

if size(sig1,1)~=size(sig2,1)
    error('sig1 and sig2 must be same length')
end
if size(sig1,1)==1
    error('sig1 and sig2 must be column vectors')
end

maxlag = (b_len-s_len)/2; % bins
segsize = maxlag*2+1; % # bins in one xcorr segment
binstart = b_len/2 - s_len/2; % zero 2nd sig before here
binend = b_len/2 + s_len/2; % zero 2nd sig after here

totalbins = size(sig1,1);
nstep_sizes = floor(totalbins/step_size - b_len/step_size) + 1;
xcos = zeros(nstep_sizes,segsize); % store each segment's xcorr
for s = 1:nstep_sizes
    starti = (s-1)*step_size + 1; % start index;
    endi = (s-1)*step_size + b_len;
    sigi = sig1(starti:endi) - mean(sig1(starti:endi));
    sigj = sig2(starti:endi) - mean(sig2(starti:endi)); % test: size should be b_len
    sigj_zeroed = sigj;
    sigj_zeroed([1:binstart-1 1+binend:b_len]) = 0; % test: sum ~=0 should be s_len)
    xco = xcorr(sigi,sigj_zeroed,maxlag);
    xcos(s,:) = xco;
end
xco_avg = mean(xcos,1);

if to_plot
    bin = 1/FS; % secs in bin, might not need
    figure
    plot(-maxlag*bin:bin:maxlag*bin,xco_avg,'k');
    xlabel('Time (sec)')
    ylabel('XCorr')
end
end