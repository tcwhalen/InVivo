function [burst] = surpriseBurst(ts, T, op)
% Tim Whalen, last edited Nov 2020
% Modification of code by Thomas Wichmann July 2005
%
% Burst detection algorithm by Legendy and Salcman (J. Neurophysiology 53:926) on
% the input ISI data stream. As the initial criterion for burst detection, spikes
% have to occur at a frequency which is higher than the baseline frequency by a factor fac.  The
% user can modify the minimal length of a prospective burst (min_length_of_burst), and the surprise
% cutoff value.  In the original paper, this value was set to 10 to detect significant spikes.
% The burst discharge is compared to a segment of data immediately preceding the burst.  A local_length
% value of 0 indicates that the entire data stream should be used to yield the firing rate (can be used
% in case of homogeneous firing of units).  Other values indicate a local time interval that should be
% used (for instance, a value of 1 would mean that each spike/burst is compared to 1 second of data
% prior to the index spike).
%
% The function uses a subroutine called surprise_new3.  The function returns 100 for
% a very high surprise value and 0 for a very low one.  This means that the output of this rou-
% tine is not accurate for such very high or very low surprise values (although bursts are correctly
% detected).
%
% The function outputs "burst" struct which contains fields describing the bursts, including the
% onset and lengths of bursts (described as indices of the input ISI stream) the highest rate within
% the burst, the average discharge rate in the burst, the pre-burst baseline rate, as well as the
% surprise values for individual bursts.  In field 1 of the structure, summary parameters are defined,
% including, num_bursts (the total number of bursts detected), mean_spikes_per_burst, total_spikes_in_bursts,
% mean_intra_burst_frequency, proportion_time_in_bursts, and proportion_spikes_in_bursts.
%
% See arguments block for input. "op" are given as name-value pairs

arguments
    %% requires
    ts (:,1) double % spike times in seconds
    T double % length of recording in seconds
    
    %% tunable parameters
    op.fac double = 1.5 % firing rate threshold (as multiple of baseline) to
        % initially consider burst
    op.min_length_of_burst double {mustBeInteger} = 2 % minimum number of 
        % ISI's a burst can have
    op.local_length double = 0 % time (in sec) in past over which to define 
        % current baseline firing rate. If 0, calculates over entire train 
        % Zero only recommended if ts is stationary; otherwise, a good 
        % choice highly depends on the type of nonstationarity
    op.surprise_cutoff double = 10 % minimum Poisson surprise for a burst to 
        % be kept
end

fac = op.fac; min_length_of_burst = op.min_length_of_burst;
local_length = op.local_length; surprise_cutoff = op.surprise_cutoff;

% vector outputs of length nbursts
burst = struct();
burst.begin = zeros(0,1);
burst.num_spikes = zeros(0,1);
burst.surprise = zeros(0,1);
burst.rate = zeros(0,1);
burst.max_rate = zeros(0,1);
burst.baseline_rate = zeros(0,1);
burst.burst_duration = zeros(0,1);
                
ISI = diff(ts);
burst_num = 0;
CA = cumsum(ISI);

if local_length == 0
    mean_FR = length(ISI)/(sum(ISI));
    fr_thr = 1/(fac*mean_FR); % calculation of frequency threshold
    beg_idx = 0;
else
    beg_idx = find(CA < local_length,1,'LAST'); % finds last index within the 'local length' - incremented by one, this will result in the first index that can be evaluate.
end
n = beg_idx;

% ***** Main loop ****

while n < length(ISI) - min_length_of_burst
    n = n+1;
    if local_length > 0
        I = ISI(find(CA > CA(n)-local_length,1,'FIRST'):n-1); % find the ISI data segment I that is fully contained within the local_length
        mean_FR = length(I)/(sum(I));
        fr_thr = 1/(fac*mean_FR); % calculation of frequency threshold
    end
    
    % ****** 1. selection step - find runs of short ISIs *******************
    if (ISI(n) < fr_thr) % finding areas of the spike train which fulfill the length_of_burst criterion
        q = 0;           % running parameter that points to the number of spikes to be added
        while (n+q <= length(ISI)) && (ISI(n+q) < fr_thr)
            q = q+1;
        end
        q = q-1;
        % at this point, the provisional burst starts at n and ends at n+q;
        % it has q+1 spikes in it
        
        % ******* 2. selection step - adjust length of burst to maximize surprise value ***************
        if q+1 >= min_length_of_burst
            m = min_length_of_burst; % try one spike more than necessary
            % running parameter of the number of spikes to be added to n
            while ((n+m <= length(ISI)) && ...
                    (ISI(n+m) < fr_thr) && ...
                    (surprise_new3(mean_FR,ISI(n:n+m)) >= surprise_new3(mean_FR,ISI(n:n+m-1)))) % 2. burst criterion - surprise values are increasing when spikes are added
                m = m+1;
            end
            m = m-1;
            
            % ******* 3. selection step - test whether adding up to 10 more ISIs will enhance surprise value **********************
            if n+m+10 <= length(ISI) % mmax is set to 10 unless one couldn't add 10 to n before reaching the end of FR
                mmax = 10;
            else
                mmax = length(ISI)-(n+m);
            end
            
            ind_long_ISI = find(ISI(n+m:n+m+mmax) > fr_thr,1,'FIRST'); % looking for 'slow spikes' within the next 10 spikes after the current burst end
            if ~isempty(ind_long_ISI) % pmax is set to be the index of the slow spike
                pmax = ind_long_ISI-1;
            else
                pmax = mmax;
            end
            
            S = surprise_new3(mean_FR,ISI(n:n+m-1));
            for p = 1:pmax % forms array of surprise values for this burst, starting from the end of the burst to pmax values later
                S(p+1) = surprise_new3(mean_FR,ISI(n:n+m-1+p));
            end
            
            if n+m < length(ISI)
                [max_S,ind_max_S] = max(S);
                if (max_S > surprise_new3(mean_FR,ISI(n:n+m-1))) % check whether the maximal surprise value in the post-burst period is higher than the maximal burst surprise
                    m = m+ind_max_S-1;
                end
            else
                m = length(ISI)-n;
            end
            
            if n+m > length(ISI)
                m = length(ISI)-n;
            end
            % at this point, the end of the index of the end of the burst
            % is settled to be n+m-1
            
            % ******** 4. selection step - test whether adjusting the front end of the burst enhances the surprise value ******************
            S = surprise_new3(mean_FR,ISI(n:n+m-1)); % surprise of unshortened train
            for p = 1: m+1-min_length_of_burst % try all possible bursts, with smallest being min_length_of_burst spikes all at end
                S(p+1) = surprise_new3(mean_FR,ISI(n+p:n+m-1));
            end
            [max_S,ind_max_S] = max(S);
            n = n+ind_max_S-1; % move start up
            m = m-ind_max_S+1; % adjust length accordingly
            
            % at this point, the beginning of the burst is settled to be n, and the length is m+1
            if (m+1 >= min_length_of_burst) && max_S > surprise_cutoff
                burst_num = burst_num + 1;
                burst.begin(burst_num) = n;
                burst.num_spikes(burst_num) = m+1;
                burst.surprise(burst_num) = max_S; 
                burst.rate(burst_num) = length(ISI(n:n+m-1))/(sum(ISI(n:n+m-1))); 
                burst.max_rate(burst_num) = 1/min(ISI(n:n+m-1));
                burst.baseline_rate(burst_num) = mean_FR;
                burst.burst_duration(burst_num) = sum(ISI(n:n+m-1));
            end
            n = n+m;
        end
    end
end

% single value outputs
burst.num_bursts = length(burst.begin);
burst.bursts_per_second = burst.num_bursts/T;
burst.mean_spikes_per_burst = mean(burst.num_spikes);
burst.median_spikes_per_burst = median(burst.num_spikes);
burst.total_spikes_in_bursts = sum(burst.num_spikes);
burst.mean_intra_burst_frequency = mean(burst.rate);
burst.median_intra_burst_frequency = median(burst.rate);
burst.proportion_time_in_bursts = (burst.total_spikes_in_bursts/burst.mean_intra_burst_frequency)/(sum(ISI(beg_idx+1:length(ISI))));
burst.proportion_spikes_in_bursts = burst.total_spikes_in_bursts/length(ISI(beg_idx+1:length(ISI)));

    function [burst,deceleration] = surprise_new3(r,data)
        % calculates surprise index
        % r = comparison firing rate (spikes per second)
        % data = ISI data to be included in the burst
        
        T2 = sum(data);
        num_spikes = length(data);
        p2 = poisscdf(num_spikes,r*T2);
        
        switch p2
            case 0 % for very small p, a value of 10exp-100 is assumed
                burst = 0;
                if nargout > 1,deceleration = 100;end
            case 1 % for very high p, a value of 1-10exp-100 is assumed
                burst = 100;
                if nargout > 1,deceleration = 0;end
            otherwise
                burst = -log(1-p2);
                if nargout > 1,deceleration = -log10(p2);end
        end
    end
end