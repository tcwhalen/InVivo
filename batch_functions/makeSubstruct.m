function [ subdata ] = makeSubstruct( data, fnum, varargin )
% Tim Whalen July 2017
% Takes data struct (see open_data) and a vector of integers such that 
% max(fnum) <= data.nfiles and returns a data struct as if only those files
% had been loaded.
% If a field of unexpected size (scalar or nfiles x 1, will display a
% warning (e.g., byunit, or fields with text information) Include an
% additional argument of 0 to turn off this printing
global verbose 
verbose = 1;
if ~isempty(varargin)
    verbose = varargin(1);
end
subdata = substructInt(data,fnum,data.nfiles);
end

function [ data2 ] = substructInt( data, fnum, nfiles )
fun = @(d) substructFun(d,fnum,nfiles);
data2 = structfun(fun,data,'UniformOutput',0);
data2.nfiles = length(fnum);
end

function [ x2 ] = substructFun( x, fnum, nfiles )
global verbose
if size(x,1) == nfiles
    x2 = x(fnum,:);
elseif size(x,2) == nfiles
    x2 = x(:,fnum);
else % if not nfiles length vector
    if isstruct(x)
        x2 = substructInt(x, fnum, nfiles); % reursion
    else
        if length(x) ~= 1 && verbose % if also not scalar
            disp('Warning: strut vector or cell array of unexpected length unchanged')
            x2 = x;
        else
            x2 = x;
        end
    end
end
end