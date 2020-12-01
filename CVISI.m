function [CV] = CVISI(ts)
% CVISI - coefficient of variation of interspike intervals
% Input: ts, vector of spike times
% Output: CV, CV of ISI's
ISI = diff(ts);
CV = std(ISI)/mean(ISI);
end

