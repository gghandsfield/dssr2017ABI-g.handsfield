function [ idx ] = gh_rangesearch( data1,data2,distance )
%gh_rangesearch Workaround to use rangesearch without Stats toolbox
%   Uses John D'Errico's ipdm function and manipulates answer to return the
%   same format as rangesearch.
%   See MATLAB documentation for function rangesearch
%
% 	Notes on usage syntax:
% 	for every row in data2, you get a column array
% 	column array tells you which rows in data1 are within distance
% 	all column arrays are stored in a cell output

distmat = ipdm(data1,data2); %distance matrix

[~,n] = size(distmat);
idx = cell(n,1);
for i = 1:n
    id = find(distmat(:,i) < distance);
    idx{i} = id;
end

end

