function [ fibertracts, elim, m, CFDfiblength ] = CFDtractography_2( XYZ, UVW, XoYoZo, STL, proxApo, distApo)
%CFDtractography This is a code to compute fibertracts from CFD data
%   XYZ is an nx3 array of coordinates, UVW is an nx3 array of vectors at
%   those coordinates, XoYoZo is an mx3 array of starting locations on the
%   surface of your mesh. STL is the STL data of your mesh. 
%   proxApo and distApo are nx3 arrays of coordinate points defining the
%   proimal and distal aponeurosis surface.
%   elim is an index of XoYoZo that did not produce fiber tracts from
%   origin to insertion.

% [IDX,~] = YiCao_knnsearch(XYZ,XoYoZo); %index of pts close to starting pts
% %remove duplicates
% [~,IA]=unique(IDX,'first');
% duplicates = setdiff(1:length(IDX),IA)';%indices of IDX that are duplicates

% better way would be to eliminate all but the closest duplicated ones,
% rather than all but the first duplicated one...
[IDX,Ds] = YiCao_knnsearch(XYZ,XoYoZo); %index of pts close to starting pts
[~,I] = sort(Ds); % sort by distance from each point
[~,I_rev] = sort(I); % reverse sort, need this for later
[~,IA] = unique(IDX(I),'first'); % indices of unique indices
dup_inds = setdiff(1:length(IDX(I)),IA)';%indices of IDX(I) that are duplicates
%temp is a logical that is 1 if it's unique and 0 if it's duplicate
temp = IDX(I); temp(IA) = 1; temp(dup_inds) = 0; 
temp = temp(I_rev); %unsort temp so it's in the same order as IDX
duplicates = find(~temp); % indices of IDX that are duplicates

m = length(IDX);
UVW = UVW./repmat(sqrt(sum(UVW.^2,2)),1,3); %unitize vector

CFDfiblength = zeros(m,1);
for i = 1:m
    
    if ismember(i,duplicates) % if the current tract is a duplicate, set length to 0 and continue
        tractpath = [];
        CFDfiblength(i) = 0;
        
    else
        [~,D1] = YiCao_knnsearch(proxApo,XYZ(IDX(i),:));
        
        if D1 > 2; % if starting point is too far from proximal Apo set length to 0 and continue
            tractpath = [];
            CFDfiblength(i) = 0;  
        else
            startpt = XYZ(IDX(i),:); tractpath = startpt;
            startdir = UVW(IDX(i),:);
            iter = 0;
            %         flip search direction if 2nd point is out of STL
            if ~inpolyhedron(STL,startpt + startdir)
                UVW = -UVW;
            end
            while inpolyhedron(STL,startpt+startdir)
                nextstep = startpt + 0.25*startdir; % 0.25mm step in startdir direction
                [newid,~] = YiCao_knnsearch(XYZ,nextstep);
                tractpath = [tractpath; nextstep];
                startpt = nextstep;
                startdir = UVW(newid,:);
                iter = iter + 1;
            end
            %         [~,D1] = YiCao_knnsearch(STL.vertices,tractpath(1,:));
            %         [~,D2] = YiCao_knnsearch(STL.vertices,tractpath(end,:));
                      
            % distance between Apo and the end of the tractpath
            [~,D2] = YiCao_knnsearch(distApo,tractpath(end,:));
            % if the distance from insertion to apo is more than 2mm, get rid of tract
            if D2 > 2
                tractpath = [];
                CFDfiblength(i) = 0;
            else
                % unsmoothed length
            %         CFDfiblength(i) = iter + D1 + D2;
            
            % smoothed length
            CFDfiblength(i) = sum(sqrt(sum((diff([smooth(tractpath(:,1))...
                smooth(tractpath(:,2)) smooth(tractpath(:,3))])).^2,2)));
            
            % linear length
            %         CFDfiblength(i) = sqrt(sum((tractpath(1,:) - tractpath(end,:)).^2,2));
            end
            
        end
        
        
    end
    
    
    
    fibertracts{i} = tractpath;
    clear tractpath
    
    
end

elim = find(CFDfiblength<10)'; %index of tracts less than 1cm
% (also includes tracts that were duplicates)

%eliminate these short tracts
fibertracts(elim) = []; 
CFDfiblength(elim) = [];

end

