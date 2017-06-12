function [ deepPennAngle, supPennAngle ] = compPennation( STL, origins, insertions, mean_vector )
%DTIpennation Compute pennation angle for DTI muscle data. Can also be used
%for CFD muscle fibers.
%   STL is the STL input data, origins is the x,y,z locations of origin
%   tracts. insertions is the x,y,z locations of insertion tracts. origin
%   and insertion should be the same length. mean_vector is the average
%   unit vector of the tract, should be the same length as origins.


% fibvecs = origins - insertions;

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
%     idx_5mm = rangesearch(facecntr,insertions,5);
    idx_5mm = gh_rangesearch(facecntr,insertions,5);

    deepPennAngle = zeros(length(idx_5mm),1);
    for rr = 1:length(idx_5mm)
        
        % mean of face normals within 5mm of each insertion
        MeanNormal = mean(STL.normals(idx_5mm{rr},:),1);
        MeanNormal = MeanNormal./norm(MeanNormal); %unitize face normal
        
        % unit vector of DTI fiber vector
%         unitvec = fibvecs(rr,:)./norm(fibvecs(rr,:));
        unitvec = mean_vector(rr,:); % mean vector of entire tract
               
        % dot product and cosine rule to determine pennation angle
        dp = dot(unitvec,MeanNormal);
        deepPennAngle(rr) = acosd( sign(dp).*dp);
    end
    deepPennAngle = abs(real(deepPennAngle) - 90);
    
%Superficial Pennation Angle (origin)
    idx_5mm = gh_rangesearch(facecntr,origins,5);
    
    % not sure what the code below is doing here. Is it needed?
%     if find(cellfun(@isempty,idx_5mm))
%         temp = find(cellfun(@isempty,idx_5mm));
%         idx_5mm(temp,1) = gh_rangesearch(facecntr,origins(temp,:),7);
%     end
    supPennAngle = zeros(length(idx_5mm),1);
    for rr = 1:length(idx_5mm)
        
        % mean of face normals within 5mm of each insertion
        MeanNormal = mean(STL.normals(idx_5mm{rr},:),1);
        MeanNormal = MeanNormal./norm(MeanNormal); %unitize face normal
        
        % unit vector of DTI fiber vector
%         unitvec = fibvecs(rr,:)./norm(fibvecs(rr,:));
        unitvec = mean_vector(rr,:); % mean vector of entire tract
        
        % dot product and cosine rule to determine pennation angle
        dp = dot(unitvec,MeanNormal);
        supPennAngle(rr) = acosd( sign(dp).*dp);
    end
    supPennAngle = abs(real(supPennAngle) - 90);

end

