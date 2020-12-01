function [ pois ] = genPoisProc( rate, T, refrac )
% Tim Whalen, last edited Nov 2020
% Generates poisson process (with optional refractory period)
% Input:
% rate: lambda for Poisson prcoess
% T: length of time to generate
% refrac: absolute refractory period (0 for genuine Poisson process)
% Output: (:,1) vector of event times

mu = 1./rate;
preal = ceil(T*rate*2); % number of exp variables to get for each unit

isis = exprnd(mu,preal,1) + refrac;
ts = cumsum(isis);
pois = ts(ts<T);
end