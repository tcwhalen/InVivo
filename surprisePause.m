function [pause] = surprisePause(ts, T, op)
% Tim Whalen, last edited Nov 2020
% Modification of code by Thomas Wichmann July 2005
%
% Am inversion of the surprise burst algorithm (Legendy and Salcman 1985)
% to instead detect "surprisingly" long sequences of ISIs (i.e. pauses).
% If you would only call something a pause if it contains no spikes at all,
% this might be more accurately called a "downstate" detector.
% A word of warning - I never found a reason to use this so it hasn't been 
% tested extensively, but pauses are an underappreciated feature of spike
% trains, so I hope you find some use for it!
%
% See arguments block for input

arguments
    %% requires
    ts (:,1) double % spike times in seconds
    T double % length of recording in seconds
    
    %% tunable parameters
    fac double = 0.3 % firing rate threshold (as multiple of baseline) to
        % initially consider pause. Unlike surpriseBurst, looking for
        % sequences UNDER this cutoff
    min_length_of_pause double {mustBeInteger} = 1 % minimum number of 
        % ISI's a pause can have (unlikely you want to change this, unless
        % you want to force pauses to contian at least some spikes)
    local_length double = 0 % time (in sec) in past over which to define 
        % current baseline firing rate. If 0, calculates over entire train 
        % Zero only recommended if ts is stationary; otherwise, a good 
        % choice highly depends on the type of nonstationarity
    surprise_cutoff double = 10 % maximum (not minimum!) Poisson surprise 
        % for a pause to be kept
end

fac = op.fac; min_length_of_pause = op.min_length_of_pause;
local_length = op.local_length; surprise_cutoff = op.surprise_cutoff;

% outputs of length npauses
pause = struct();
pause.begin = zeros(0,1);
pause.num_spikes = zeros(0,1);
pause.surprise = zeros(0,1);
pause.rate = zeros(0,1);
pause.max_rate = zeros(0,1);
pause.baseline_rate = zeros(0,1);
pause.pause_duration = zeros(0,1);
                
ISI = diff(ts);
pause_num = 0;
CA = cumsum(ISI);

if local_length == 0
    mean_FR = length(ISI)/(sum(ISI));
    fr_thr = 1/(fac*mean_FR);   % calculation of frequency threshold
    beg_idx = 0;
else
    beg_idx = find(CA < local_length,1,'LAST');     % finds last index within the 'local length' - incremented by one, this will result in the first index that can be evaluate.
end
n = beg_idx;

% ***** Main loop ****

while n < length(ISI) - min_length_of_pause
    n = n+1;
    if local_length > 0
        I = ISI(find(CA > CA(n)-local_length,1,'FIRST'):n-1);     % find the ISI data segment I that is fully contained within the local_length
        mean_FR = length(I)/(sum(I));
        fr_thr = 1/(fac*mean_FR);                                   % calculation of frequency threshold
    end
    
    % ****** 1. selection step - find runs of long ISIs *******************
    if (ISI(n) > fr_thr) % finding areas of the spike train which fulfill the length_of_pause criterion
        q = 0;           % running parameter that points to the number of spikes to be added
        while (n+q <= length(ISI)) && (ISI(n+q) > fr_thr)
            q = q+1;
        end
        q = q-1;
        % at this point, the provisional pause starts at n and ends at n+q;
        % it has q+1 spikes in it
        
        % ******* 2. selection step - adjust length of pause to minimize surprise value ***************
        if q+1 >= min_length_of_pause
            
            m = min_length_of_pause; % try one spike more than necessary
            % running parameter of the number of spikes to be added to n
            while ((n+m <= length(ISI)) && ...
                    (ISI(n+m) > fr_thr) && ...
                    (surprise_new3(mean_FR,ISI(n:n+m)) <= surprise_new3(mean_FR,ISI(n:n+m-1)))),   % 2. pause criterion - surprise values are increasing when spikes are added
                m = m+1;
            end
            m = m-1;
            % at this point, the beginning of the pause is provisionally settled to be n, the end n+m
            % the pause has m+1 spikes in it.
            
            % ******* 3. selection step - test whether adding up to 10 more ISIs will diminish surprise value **********************
            if n+m+10 <= length(ISI) % mmax is set to 10 unless one couldn't add 10 to n before reaching the end of FR
                mmax = 10;
            else
                mmax = length(ISI)-(n+m);
            end
            
            ind_long_ISI = find(ISI(n+m:n+m+mmax) < fr_thr,1,'FIRST'); % looking for 'fast spikes' within the next 10 spikes after the current pause end
            if ~isempty(ind_long_ISI) % pmax is set to be the index of the slow spike
                pmax = ind_long_ISI-1;
            else
                pmax = mmax;
            end
            
            S = surprise_new3(mean_FR,ISI(n:n+m-1));
            for p = 1:pmax % forms array of surprise values for this pause, starting from the end of the pause to pmax values later
                S(p+1) = surprise_new3(mean_FR,ISI(n:n+m-1+p));
            end
            
            if n+m < length(ISI)
                [min_S,ind_min_S] = min(S);
                if (min_S < surprise_new3(mean_FR,ISI(n:n+m-1))) % check whether the min surprise value in the post-pause period is lower than the min pause surprise
                    m = m+ind_min_S-1;
                end
            else
                m = length(ISI)-n; % don't think this can ever be reached?
            end
            
            if n+m > length(ISI)
                m = length(ISI)-n; % also seems impossible to reach
            end
            % at this point, the end of the index of the end of the pause
            % is settled to be n+m-1
            
            % ******** 4. selection step - test whether adjusting the front end of the pause enhances the surprise value ******************
            S = surprise_new3(mean_FR,ISI(n:n+m-1)); % surprise of unshortened train
            for p = 1: m+1-min_length_of_pause % try all possible pauses, with smallest being min_length_of_pause spikes all at end
                S(p+1) = surprise_new3(mean_FR,ISI(n+p:n+m-1));
            end
            [min_S,ind_min_S] = min(S);
            n = n+ind_min_S-1; % move start up
            m = m-ind_min_S+1; % adjust length accordingly
            
            % at this point, the beginning of the pause is settled to be n, and the length is m+1
            if (m+1 >= min_length_of_pause) && min_S < surprise_cutoff 
                pause_num = pause_num + 1;
                pause.begin(pause_num) = n;
                pause.num_spikes(pause_num) = m+1;
                pause.surprise(pause_num) = min_S; % TW: again no need to recalculate
                pause.rate(pause_num) = length(ISI(n:n+m-1))/(sum(ISI(n:n+m-1))); % TW: m -> m-1
                pause.max_rate(pause_num) = 1/min(ISI(n:n+m-1)); % TW: m -> m-1
                pause.baseline_rate(pause_num) = mean_FR;
                pause.pause_duration(pause_num) = sum(ISI(n:n+m-1));
            end
            n = n+m;
        end
    end
end

pause.num_pauses = length(pause.begin);
pause.pauses_per_second = pause.num_pauses/T;
pause.mean_spikes_per_pause = mean(pause.num_spikes);
pause.median_spikes_per_pause = median(pause.num_spikes);
pause.total_spikes_in_pauses = sum(pause.num_spikes);
pause.mean_intra_pause_frequency = mean(pause.rate);
pause.median_intra_pause_frequency = median(pause.rate);
pause.proportion_time_in_pauses = (pause.total_spikes_in_pauses/pause.mean_intra_pause_frequency)/(sum(ISI(beg_idx+1:length(ISI))));
pause.proportion_spikes_in_pauses = pause.total_spikes_in_pauses/length(ISI(beg_idx+1:length(ISI)));

    function [pause,deceleration] = surprise_new3(r,data)
        % calculates surprise index
        % r = comparison firing rate (spikes per second)
        % data = ISI data to be included in the pause
        
        T2 = sum(data);
        num_spikes = length(data);
        p2 = poisscdf(num_spikes,r*T2);
        
        switch p2
            case 0 % for very small p, a value of 10exp-100 is assumed
                pause = 0;
                if nargout > 1,deceleration = 100;end
            case 1 % for very high p, a value of 1-10exp-100 is assumed
                pause = 100;
                if nargout > 1,deceleration = 0;end
            otherwise
                pause = -log(1-p2);
                if nargout > 1,deceleration = -log10(p2);end
        end
    end
end