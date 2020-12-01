function [ind] = find_local_nmax(X,n)
% [ind] = local_nmax(X,n)
% This program finds all peaks in vector X
% and returns their index position
%
%	Usage:  [ind] = find_local_max(X)
%		where 	X = vector of arbitrary length
%				n = find max over n samples (must be odd)
%				ind = indices into X where peak occur
% Created by RST, 2003-11-2

sz = size(X);
long_dims = find(sz>1);
if length(sz)>2 | length(long_dims)>1
	error('find_local_nmax doesn''t know how to work on arrays');
elseif sz(2)>1
	X = X';
end
if (n/2 - fix(n/2)) == 0
	error('n must be an odd number when using find_local_nmax');
end


% First pad ends of data (by mirroring)
data_len = length(X);
filt_len = (n-1)/2;
filt_mid = filt_len+1;

if max(sz)<filt_len
	ind = [];
	return
end

filt_ind = (-filt_len:filt_len);

start_ra = X( filt_len:-1:1);
stop_ra = X( end:-1:(end-filt_len+1));
	
pad_X = [start_ra; X; stop_ra];
start = filt_len+1;			%start of real data
stop = data_len+start-1;	%stop of real data

j = 1;
for i = start:stop
	if pad_X(i) == max( pad_X(i+filt_ind) )
		ind(j) = i-filt_len;
		j = j+1;
	end
end

% Peaks should be >n points apart
n_inds = length(ind);
all = ones(1,n_inds);
good_inds = [];
for i = 1:n_inds
	
	not_i = logical(all);
	not_i(i) = false;
	near_log = abs( ind(i)-ind(not_i))<filt_len ;
	near_inds = ind(not_i);
	
	% save if no near peaks
	if ~any(near_log)
		good_inds = [good_inds ind(i)];
	else
		% save if highest of nearby peaks
		if X( ind(i) ) > max( X( near_inds(near_log) ))
			good_inds = [good_inds ind(i)];
			
		else
		
			% save index if closest to the mid-point of equal nearby peaks
			equal_log = ( X( near_inds) == X(ind(i)) )' ;
			if any( equal_log & near_log )
				mid = round(mean( [ind(i) near_inds(equal_log & near_log )] ));

				near_difs = mid - near_inds(equal_log & near_log );
				if abs(mid-ind(i)) < min( abs(near_difs) )	% ind(i) is closer to mid
					good_inds = [good_inds ind(i)];
				elseif (abs(mid-ind(i)) == min(abs(near_difs))) & (mid-ind(i))>0
					good_inds = [good_inds ind(i)];
				end
			end
		end
	end
end
ind = good_inds;


