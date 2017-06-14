function [ deepPennAngle, supPennAngle ] = compPennation( STL, origins, insertions, mean_vector )
%compPennation Compute pennation angles for muscle data. 
%   INPUTS:
%   STL is the STL muscle input data
%   origins is x,y,z locations of n origin tracts (nx3)
%   insertions is x,y,z locations of n insertion tracts (nx3)
%   mean_vector is average unit vector of n tracts (nx3)
%
%   OUTPUTS:
%   deepPennAngle is angles between tract and distal aponeurosis (nx1)
%   supPennAngle is angle between tract and prox aponeurosis (nx1)

%unitize mean vector
unitvec = mean_vector./sqrt(sum(mean_vector.^2,2));

% find center of every stl face
 facecntr = zeros(length(STL.faces),3);
    for qq = 1:length(STL.faces)
        
        %vertices surrounding each STL face
        verts = [STL.vertices(STL.faces(qq,:),1),...
            STL.vertices(STL.faces(qq,:),2),...
            STL.vertices(STL.faces(qq,:),3)];
        
        %center of every face in STL
        facecntr(qq,:) = mean(verts,1); 
    end
    
% Deep Pennation Angle (insertion)
    idx_5mm = gh_rangesearch(facecntr,insertions,5); %indices of facecenter within 5mm of insertions

    deepPennAngle = zeros(length(idx_5mm),1);
    for rr = 1:length(idx_5mm)
        
        % mean of face normals within 5mm of each insertion
        MeanNormal = mean(STL.normals(idx_5mm{rr},:),1);
        MeanNormal = MeanNormal./norm(MeanNormal); %unitize face normal
           
        % dot product and cosine rule to determine pennation angle
        dp = dot(unitvec(rr,:),MeanNormal);
        deepPennAngle(rr) = acosd( sign(dp).*dp);
    end
    deepPennAngle = abs(real(deepPennAngle) - 90);
    
%Superficial Pennation Angle (origin)
    idx_5mm = gh_rangesearch(facecntr,origins,5); %indices of facecenter within 5mm of origins
    
    % code below is to deal with origins that have no facecenters within 5mm, reassess at 7mm
%     if find(cellfun(@isempty,idx_5mm))
%         temp = find(cellfun(@isempty,idx_5mm));
%         idx_5mm(temp,1) = gh_rangesearch(facecntr,origins(temp,:),7);
%     end
    supPennAngle = zeros(length(idx_5mm),1);
    for rr = 1:length(idx_5mm)
        
        % mean of face normals within 5mm of each insertion
        MeanNormal = mean(STL.normals(idx_5mm{rr},:),1);
        MeanNormal = MeanNormal./norm(MeanNormal); %unitize face normal
              
        % dot product and cosine rule to determine pennation angle
        dp = dot(unitvec(rr,:),MeanNormal);
        supPennAngle(rr) = acosd( sign(dp).*dp);
    end
    supPennAngle = abs(real(supPennAngle) - 90);

end

