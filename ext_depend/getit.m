function c=getit(c)
% like cell2mat but cell 2 smaller cell
if iscell(c)
    c = cellfun(@getit, c, 'UniformOutput', 0);
    c = cat(2,c{:});
else
    c = {c};
end

end