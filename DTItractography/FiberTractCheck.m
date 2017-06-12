function [ idxOut, LengthOut, neworigin, newinsertion] ...
    = FiberTractCheck( DataIn, OriginVerts, InsertionVerts, STL )
%FiberTractCheck This function sorts through a set of fiber tracts and
%finds those whose origins and insertions are within 1mm of the surface
%vertices.
%   DataIn is the data structure from DSIStudio, SurfaceVerts is the
%   vertices from stlread, idxOut is the fiber index of the fibers that
%   meet the criteria, LengthOut is the lengths of those fibres. D1 is the
%   distance between the origin point and the STL. neworigin are the points
%   on the proximal aponeurosis closest to the input origin. newinsertion
%   are the points on the distal aponeurosis closest to the input insertion

data = DataIn;

origins = data.tracts_xyz(:,data.fibindex(:,1))';
insertions = data.tracts_xyz(:,data.fibindex(:,2))';

% switch origins and insertions if the origin is distal to insertion
testvec = insertions - origins;
testidx = find(testvec(:,3) > 0);
testor = origins(testidx,:);
testins = insertions(testidx,:);
origins(testidx,:) = testins;
insertions(testidx,:) = testor;

%% Sort Origin and Insertion Points by distance from [0,0,0]
% test = sqrt(sum(OriginVerts.^2,2));
% [~,I] = sort(test);
% OriginVerts = OriginVerts(I,:);
% 
% test = sqrt(sum(InsertionVerts.^2,2));
% [~,I] = sort(test);
% InsertionVerts = InsertionVerts(I,:);

%%
[or_id,D1] = YiCao_knnsearch(OriginVerts,origins); % nearest surface points to tract origins
[ins_id,D2] = YiCao_knnsearch(InsertionVerts,insertions); % nearest surface points to tract insertions
neworigin = OriginVerts(or_id,:); % points that are actually on aponeurosis
newinsertion = InsertionVerts(ins_id,:); % points that are actually on aponeurosis

%tract idx where origin and insertion are both within 1mm of aponeuroses
keepers = intersect(find(D1 < 1),find(D2 < 1)); 

%% simple way to figure out extra length between origin and insertion not captured in tract length

inout_orig = inpolyhedron(STL,origins(keepers,:)); % logical of whether original origin is inside STL
distsign_o = inout_orig + inout_orig - 1; % +1/-1 array from logical array

inout_ins = inpolyhedron(STL,insertions(keepers,:));% logical of whether original origin is inside STL
distsign_i = inout_ins + inout_ins - 1; % +1/-1 array from logical array

%min dist between origin/insertion and STL with a sign indicating whether
%it's inside or outside of the muscle
% ExtraLength = distsign_o.*D1(keepers) + distsign_i.*D2(keepers); 

%% Alternative Way: Sample points along line from tract origin to insertion, find point closest to aponeurosis
origins = origins(keepers,:);
insertions = insertions(keepers,:);
testvec = insertions - origins;
testvec = testvec./ repmat(sqrt(sum(testvec.^2,2)),1,3);

extradist_o = zeros(length(keepers),1);
extradist_i = extradist_o;
for i = 1:length(keepers)
    alongpts_o = repmat(origins(i,:),41,1) + repmat((0:0.05:2)',1,3) .* ...
        repmat(-distsign_o(i) * testvec(i,:),41,1);
    odmat = ipdm(OriginVerts,alongpts_o);
    [~,b] = min(min(odmat,[],1)); % b is the index of the close point
    neworigin(i,:) = alongpts_o(b,:); % pt that is very close to aponeurosis
    % distance between new close point and the original tract endpoint
    extradist_o(i) = sqrt(sum((alongpts_o(b,:) - alongpts_o(1,:)).^2,2));
    
    alongpts_i = repmat(insertions(i,:),41,1) + repmat((0:0.05:2)',1,3) .* ...
        repmat(distsign_o(i) * testvec(i,:),41,1);
    idmat = ipdm(InsertionVerts,alongpts_i);
    [~,b] = min(min(idmat,[],1)); % b is the index of the close point
    newinsertion(i,:) = alongpts_i(b,:); % pt that is very close to aponeurosis
    
    % distance between new close point and the original tract endpoint
    extradist_i(i) = sqrt(sum((alongpts_i(b,:) - alongpts_i(1,:)).^2,2));
end

%distance between extrapolated origin/insertion and the one on the DTI
%tract. This should be added to the tract length;
ExtraLength = distsign_o.*extradist_o + distsign_i.*extradist_i;


%% Outputs
idxOut = keepers; % fiber index of fibers that meet criteria

% % ArcLength
LengthOut = data.length(keepers)' + ExtraLength;


neworigin = neworigin(keepers,:);
newinsertion = newinsertion(keepers,:);

% linear length
% LengthOut = sqrt(sum((neworigin - newinsertion).^2,2));
    
    
end

