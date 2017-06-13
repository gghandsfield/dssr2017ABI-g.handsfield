function [ idx ] = gh_rangesearch( data1,data2,distance )
%gh_rangesearch Workaround to use rangesearch without Stats toolbox
%   Uses John D'Errico's ipdm function and manipulates answer to return the
%   same format as rangesearch.
%   See MATLAB documentation for function rangesearch

distmat = ipdm(data1,data2); %distance matrix

[~,n] = size(distmat);
idx = cell(n,1);
for i = 1:n
    id = find(distmat(:,i) < distance);
    idx{i} = id;
end

end

