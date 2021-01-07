function [ CV2, CV2_dist ] = CV2ISI( ts )
% Tim Whalen January 2021
% CV_2 spike regularity metric as described in Holt et al. 1996, "Comparison
% of discharge variability in vitro and in vivo in cat visual cortex neurons"
% which better controls for nonstationarities in firing rate
% Inputs: ts, (:,1) vector of spike times
% Outputs: CV2, the mean of the CV2 distribution
%          CV2_dist: (:,1), the full CV2 distribution (usually not needed)

ISI = diff(ts);
CV2_dist = zeros(length(ISI)-1,1);
for i = 1:length(ISI)-1
    k = ISI(i);
    j = ISI(i+1);
    CV2_dist(i) = (2*abs(k-j))/(k+j);
end
CV2 = mean(CV2_dist);