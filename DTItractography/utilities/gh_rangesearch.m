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

% % Old Method
% [a,b] = find(distmat < distance); % row and column index of locations within distance
% 
% sorted = sortrows([b a]); % sort the row and column indices by column
% steps = find(diff(sorted(:,1))); % find locations where row changes
% 
% %for each row, assemble corresponding columns into cell array
% idx = cell(length(steps),1); 
% begin = 1;
% for i = 1:length(steps)
%     idx{i,1} = sorted(begin:steps(i),2)';
%     begin = steps(i) + 1;
% end
% idx{i+1,1} = sorted(begin:end,2)';


end

