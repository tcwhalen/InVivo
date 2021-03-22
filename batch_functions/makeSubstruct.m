function [ subdata ] = makeSubstruct( data, fnum, varargin )
% Tim Whalen July 2017
% Takes data struct (see open_data) and a vector of integers such that 
% max(fnum) <= data.nfiles and returns a data struct as if only those files
% had been loaded.
% If a field is of unexpected size (i.e. not scalar or (nfiles,1), such as
% a byunit field or text information), these fields will be copied in full 
% and a warning will be printed. Include an additional argument of 0 to 
% turn off these warnings.

global verbose 
verbose = 1;
if ~isempty(varargin)
    verbose = varargin(1);
end
subdata = substructInternal(data,fnum,data.nfiles);
end

function [ data2 ] = substructInternal( data, fnum, nfiles )
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
        x2 = substructInternal(x, fnum, nfiles); % reursion
    else
        if length(x) ~= 1 && verbose % if also not scalar
            disp('Warning: strut vector or cell array of unexpected length unchanged')
        end
        x2 = x;
    end
end
end