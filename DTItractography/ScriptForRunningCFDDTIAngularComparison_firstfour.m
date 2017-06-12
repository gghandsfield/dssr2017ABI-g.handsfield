%% This script is to find the angular differences between CFD and DTI data in the first 4 subjects
% Use this script for no FlowGuides

%% Declaring the Output Variables
CFDtracts = cell(1);
DTItracts = cell(1);
cfdOut = cell(1);
dtiOut = cell(1);
FibLengths = cell(1);
FibPennation = cell(1);
prop30 = cell(1);
prop20 = cell(1);
prop10 = cell(1);
    
%% Main Script Initiation
ui = input(['Analyze Training Models with\n'...
    '1."no flow guides" or\n'...
    '2."flow guides"?\n']);

% Pointers to CFD nodal folders
    CFDpre = 'D:\Users\ghan078\Desktop\DTI-CFD_MGComp\';

    if ui == 1;
        CFDfolders = {'Male1_InvSeg_nostl_80percentProxApo_2SimpleFlowGuides_edited_PM_edit_FG_off\Design 1\Scenario 1\solver'...
            'Female2_InvSeg_80percentProxApo_2SimpleFlowGuides_edited_FG_off\Design 1\Scenario 1\solver'...
            'Male3_InvSeg_shaperedo_80percentProxApo_2SimpleFlowGuides_AreaThirds_FG_off\Design 1\Scenario 1\solver'...
            'Female4_InvSeg_80percentProxApo_2SimpleFlowGuides_edited_2_FG_off\Design 1\Scenario 1\solver'};
        myTitle = 'no flow guides';
        
    elseif ui == 2;
        CFDfolders = {'Male1_InvSeg_nostl_80percentProxApo_2SimpleFlowGuides_edited_PM_edit_FG_on\Design 1\Scenario 1\solver'...
            'Female2_InvSeg_80percentProxApo_2SimpleFlowGuides_edited_FG_on\Design 1\Scenario 1\solver'...
            'Male3_InvSeg_shaperedo_80percentProxApo_2SimpleFlowGuides_AreaThirds_FG_on_slightchange\Design 1\Scenario 1\solver'...
            'Female4_InvSeg_80percentProxApo_2SimpleFlowGuides_edited_2_FG_on_slightchange3\Design 1\Scenario 1\solver'};
        myTitle = 'flow guides';
    
    end



for ii = 1:4
    loc = ['D:\Users\ghan078\Documents\corrected_data\Subj0' num2str(ii)];
    
    cd(loc); P = pwd; subjind = P(end-5:end);
    
    %% Setup DTI coordinates and vectors
    % directory for DTI tract data
    cd(['D:\Users\ghan078\Documents\DTI_tracts\' subjind]);
    data = load('DTItracts1');
    data2 = load('DTItracts2');
    
    %DTI coordinates and vectors
    [vecs,XYZmid] = DTIvecsfromTracts(data);
    DTIdata1 = [XYZmid vecs];
    [vecs,XYZmid] = DTIvecsfromTracts(data2);
    DTIdata2 = [XYZmid vecs];
    
    %% Setup CFD coordinates and vectors
    
    loc = [CFDpre CFDfolders{ii}]; % indexed folders
    cd(loc);
    M = csvread('Scenario 1_nodal.csv',1,1);
    CFDcoord = M(:,2:4);
    CFDvect = M(:,5:7);
    % eliminate rows of all zeros
    indrow = find(CFDcoord(:,1) == 0 & CFDcoord(:,2) == 0 & CFDcoord(:,3) == 0);
    CFDcoord(indrow,:) = []; CFDvect(indrow,:) = [];
    
    % rigid transformation if necessary (b/c coordinate frame changed)
    cd(P); % change back to the original folder
    if exist('transformationdata.mat','file') %load a previously saved transformation if it exists
        load('transformationdata.mat'); %Z = b*MAT*T + c;
        %trans.c is translation; trans.T is rotation; trans.b is scaling;
        
        [STL.vertices, STL.faces, STL.normals] = stlread('MG_surf.stl'); % new STL info
    else
        [STL.vertices, STL.faces, STL.normals] = stlread('MG_surf.stl'); % new STL info
        cd(['D:\Users\ghan078\Documents\for_geoff\for_geoff\' subjind]) % switch to old directory
        [oldSTL.vertices, oldSTL.faces] = stlread('MG_surf.stl'); % old STL info
        cd(P); % switch back to original folder
        % perform procrustes and save file if it isn't on the path
        [~,~,trans] = procrustes(STL.vertices,oldSTL.vertices);
%         save('transformationdata.mat','trans');
        transvertices = trans.b*oldSTL.vertices*trans.T + trans.c; %transform old STL
        
        % plot the old STL in red and new one in black to examine overlap
        figure; hold on;
        patch('vertices',STL.vertices ,'faces',STL.faces,...
            'facealpha',0.3,'facecolor','none','edgecolor','k');
        patch('vertices',transvertices ,'faces',oldSTL.faces,...
            'facealpha',0.3,'facecolor','none','edgecolor','r');
    end
    % fix length of trans.c so it can be used with CFDcoord and CFDvect
    c = repmat(trans.c(1,:),length(CFDcoord),1);
    %Transform CFD info into new DTI coordinate frame
    CFDcoord = trans.b*CFDcoord*trans.T + c;
    CFDvect = trans.b*CFDvect*trans.T; %vectors don't get shifted, just rotated
    CFDunit = CFDvect./repmat(sqrt(sum(abs(CFDvect).^2,2)),1,3); % CFD unit vector
    
    CFDdata = [CFDcoord CFDunit];
    
    %% Vector Compare
    
    [theta, CFDunit, CFDcoord, DTIunit, DTIcoord] = VectorCompare_working(CFDdata,DTIdata1,DTIdata2);
    
    % decide whether to use data1 or data2
    if numel(find(theta{1,1}<10))/length(theta{1,1}) > numel(find(theta{1,2}<10))/length(theta{1,2})
        idx = 1;
    else
        idx = 2;
        data = data2;
    end


    %% Load Inventor Aponeuroses and Transform Coordinate Frame
    [proxSTL.v, proxSTL.f] = stlread('proxApo_Inventor.stl');
    [distSTL.v, distSTL.f] = stlread('distApo_Inventor.stl');
    
    %PROXIMAL APO
    % fix length of c
    c = repmat(trans.c(1,:),length(proxSTL.v),1);
    %Transform STL verts info into new DTI coordinate frame
    proxSTL.v = trans.b*proxSTL.v*trans.T + c;
    
    %DISTAL APO
    % fix length of c
    c = repmat(trans.c(1,:),length(distSTL.v),1);
    %Transform STL verts info into new DTI coordinate frame
    distSTL.v = trans.b*distSTL.v*trans.T + c;

    figure; hold on;
        patch('vertices',STL.vertices ,'faces',STL.faces,...
            'facealpha',0.3,'facecolor','none','edgecolor','k');
        patch('vertices',proxSTL.v ,'faces',proxSTL.f,...
            'facealpha',0.3,'facecolor','none','edgecolor','r');
        patch('vertices',distSTL.v ,'faces',distSTL.f,...
            'facealpha',0.3,'facecolor','none','edgecolor','r');
    
    %% DTI Muscle Architecture
    
    % index and fiber length of tracts that run from origin to insertion
    [idx1, DTIfiblength, neworigin, newinsertion] = FiberTractCheck(data, proxSTL.v, distSTL.v, STL);    
        
    % start and end point of all tracts from DSIStudio
    allOrigins = data.tracts_xyz(:,data.fibindex(:,2))';
    allInsertions = data.tracts_xyz(:,data.fibindex(:,1))';
    
    % origins and insertions that pass tests in FiberTractCheck
    origins = allOrigins(idx1,:);
    insertions = allInsertions(idx1,:);
    
    % mean vectors of each tract
    meanvec_dti = zeros(length(idx1),3);
    for w = 1:length(idx1)
        tract = data.tracts_xyz(:,data.fibindex(idx1(w),1):data.fibindex(idx1(w),2))';
        vec = diff(tract);
        meanvec_dti(w,:) = mean(vec);
    end
    
    % Pennation from DTI
    [deepDTIpenn, supDTIpenn] = compPennation(STL,neworigin,newinsertion, meanvec_dti);
    
    %% CFD
    %     load('orgncoords.mat'); % a set of manually selected pts on prox apo
    % use 'origins' which is the origins pulled from DTI tractography
    [fibertracts_CFD, elim, m, CFDfiblength] = CFDtractography( CFDcoord{idx}, CFDunit{idx}, origins, STL, proxSTL.v, distSTL.v);
    
    cfd_origins = zeros(1,3);
    cfd_insertions = cfd_origins;
    meanvec_cfd = zeros(length(fibertracts_CFD),3);
    for w = 1:length(fibertracts_CFD)
        cfd_origins(w,:) = fibertracts_CFD{w}(1,:);
        cfd_insertions(w,:) = fibertracts_CFD{w}(end,:);
        meanvec_cfd(w,:) = mean(diff(fibertracts_CFD{w}));
    end
    meanvec_cfd = meanvec_cfd./ repmat(sqrt(sum(meanvec_cfd.^2,2)),1,3);
    
    %Pennation from CFD
    [deepCFDpenn, supCFDpenn] = compPennation(STL,cfd_origins,cfd_insertions, meanvec_cfd);
    
    %% Comparing DTI and CFD
    
    %remove DTI  fibers that have no corresponding fiber in CFD tract set
    
    idx1(elim) = [];
    DTIfiblength(elim) = [];
    deepDTIpenn(elim) = [];
    supDTIpenn(elim) = [];
    origins(elim,:) = [];
    insertions(elim,:) = [];
    
    %Plot CFD tracts
    figure; hold on
    for j = 1:length(fibertracts_CFD)
%         plot3(fibertracts_CFD{j}(:,1),fibertracts_CFD{j}(:,2),fibertracts_CFD{j}(:,3))
        plot3(smooth(fibertracts_CFD{j}(:,1)),smooth(fibertracts_CFD{j}(:,2)),...
            smooth(fibertracts_CFD{j}(:,3)),'-k');
    end
% %        linear
%     for j = 1:length(fibertracts_CFD)
%         %         plot3(fibertracts_CFD{j}(:,1),fibertracts_CFD{j}(:,2),fibertracts_CFD{j}(:,3))
%         plot3([fibertracts_CFD{j}(1,1);fibertracts_CFD{j}(end,1)],...
%             [fibertracts_CFD{j}(1,2);fibertracts_CFD{j}(end,2)],...
%             [fibertracts_CFD{j}(1,3);fibertracts_CFD{j}(end,3)]);
%     end
    axis equal
    title(['Subj0' num2str(ii) 'CFD tracts' myTitle])
    patch('vertices',STL.vertices ,'faces',STL.faces,...
        'LineStyle','none',...
        'facealpha',0.1,'facecolor','r','edgecolor','k');
%     savefig(['Subj0' num2str(ii) 'CFD tracts' myTitle '.fig'])
    
    % Plot linear DTI tracts
    %lineartracts is origins interleaved with insertions
    lineartracts = reshape([origins(:) insertions(:)]',2*size(origins,1),[]);
    figure; hold on
    for j = 1:length(origins)
        plot3(lineartracts(2*j-1:2*j,1),lineartracts(2*j-1:2*j,2),lineartracts(2*j-1:2*j,3));
    end
    axis equal
    title(['Subj0' num2str(ii) 'DTI linear tracts'])
    patch('vertices',STL.vertices ,'faces',STL.faces,...
        'LineStyle','none',...
        'facealpha',0.1,'facecolor','r','edgecolor','k');
%     savefig(['Subj0' num2str(ii) 'DTI linear tracts.fig'])
    
    % Plot DTI tracts
    figure; hold on
    fibertracts_DTI = cell(1);
    for j = 1:length(origins);
        plot3(data.tracts_xyz(1,data.fibindex(idx1(j),1):data.fibindex(idx1(j),2)),...
            data.tracts_xyz(2,data.fibindex(idx1(j),1):data.fibindex(idx1(j),2)),...
            data.tracts_xyz(3,data.fibindex(idx1(j),1):data.fibindex(idx1(j),2)),'-k');
        fibertracts_DTI{j} = data.tracts_xyz(:,data.fibindex(idx1(j),1):data.fibindex(idx1(j),2));
    end
    axis equal; title(['Subj0' num2str(ii) 'DTI tracts'])
    patch('vertices',STL.vertices ,'faces',STL.faces,...
        'LineStyle','none',...
        'facealpha',0.1,'facecolor','r','edgecolor','k');
%     savefig(['Subj0' num2str(ii) 'DTI tracts' myTitle '.fig'])
    
    
    %% Setting up some plots
    
    reds = find(theta{1,idx}>30); % locations where theta is greater than 30 degrees
    blues = find(theta{1,idx}<=10); % locations where theta is less than 10 degrees
    cyans = find(theta{1,idx}<=20 & theta{1,idx}>10);% locations where theta is less than 20 degrees
    greens = find(theta{1,idx}<=30 & theta{1,idx}>20);% locations where theta is less than 30 degrees
    
    % this makes nx3 array that is non-zero only for the vectors that
    % correspond with the accuracy of the specified color
    CFDzeros = zeros(length(CFDunit{idx}),3);
    CFDreds = CFDzeros; CFDreds(reds,1:3) = CFDunit{idx}(reds,:);
    CFDblues = CFDzeros; CFDblues(blues,1:3) = CFDunit{idx}(blues,:);
    CFDcyans = CFDzeros; CFDcyans(cyans,1:3) = CFDunit{idx}(cyans,:);
    CFDgreens = CFDzeros; CFDgreens(greens,1:3) = CFDunit{idx}(greens,:);
    
    % this makes nx3 array that is non-zero only for the vectors that
    % correspond with the accuracy of the specified color
    DTIzeros = zeros(length(DTIunit{idx}),3);
    DTIreds = DTIzeros; DTIreds(reds,1:3) = DTIunit{idx}(reds,:);
    DTIblues = DTIzeros; DTIblues(blues,1:3) = DTIunit{idx}(blues,:);
    DTIcyans = DTIzeros; DTIcyans(cyans,1:3) = DTIunit{idx}(cyans,:);
    DTIgreens = DTIzeros; DTIgreens(greens,1:3) = DTIunit{idx}(greens,:);
    
    %% Plotting color coded CFD vectors
    figure; hold on;
    quiver3(CFDcoord{idx}(:,1),CFDcoord{idx}(:,2),CFDcoord{idx}(:,3),...
        CFDreds(:,1),CFDreds(:,2),CFDreds(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[1 0 0]);
    quiver3(CFDcoord{idx}(:,1),CFDcoord{idx}(:,2),CFDcoord{idx}(:,3),...
        CFDblues(:,1),CFDblues(:,2),CFDblues(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[0 0 1]);
    quiver3(CFDcoord{idx}(:,1),CFDcoord{idx}(:,2),CFDcoord{idx}(:,3),...
        CFDcyans(:,1),CFDcyans(:,2),CFDcyans(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[0 1 1]);
    quiver3(CFDcoord{idx}(:,1),CFDcoord{idx}(:,2),CFDcoord{idx}(:,3),...
        CFDgreens(:,1),CFDgreens(:,2),CFDgreens(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[0 1 0]);
    axis equal
    title(['Subj0' num2str(ii) 'arrows' myTitle]);
%     savefig(['Subj0' num2str(ii) 'CFDarrows_' myTitle '.fig']);
    
    %% Plotting color coded DTI vectors
    figure; hold on;
    quiver3(DTIcoord{idx}(:,1),DTIcoord{idx}(:,2),DTIcoord{idx}(:,3),...
        DTIreds(:,1),DTIreds(:,2),DTIreds(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[1 0 0]);
    quiver3(DTIcoord{idx}(:,1),DTIcoord{idx}(:,2),DTIcoord{idx}(:,3),...
        DTIblues(:,1),DTIblues(:,2),DTIblues(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[0 0 1]);
    quiver3(DTIcoord{idx}(:,1),DTIcoord{idx}(:,2),DTIcoord{idx}(:,3),...
        DTIcyans(:,1),DTIcyans(:,2),DTIcyans(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[0 1 1]);
    quiver3(DTIcoord{idx}(:,1),DTIcoord{idx}(:,2),DTIcoord{idx}(:,3),...
        DTIgreens(:,1),DTIgreens(:,2),DTIgreens(:,3),'AutoScale','on','AutoScaleFactor',1,'Color',[0 1 0]);
    axis equal
    title(['Subj0' num2str(ii) 'tracts' num2str(idx) 'DTIarrows_' myTitle]);
%     savefig(['Subj0' num2str(ii) 'DTIarrows_' myTitle '.fig']);
    
    
    p30 = (numel(blues)+numel(cyans)+numel(greens))/(numel(blues)+numel(cyans)+numel(greens)+numel(reds));
    p20 = (numel(blues)+numel(cyans))/(numel(blues)+numel(cyans)+numel(greens)+numel(reds));
    p10 = (numel(blues))/(numel(blues)+numel(cyans)+numel(greens)+numel(reds));
    
    %% Some Output Variables
    CFDtracts{ii,1} = fibertracts_CFD';
    DTItracts{ii,1} = fibertracts_DTI';
    
    cfdOut{ii} = [CFDcoord{idx} CFDunit{idx}];
    dtiOut{ii} = [DTIcoord{idx} DTIunit{idx}];
    
    FibLengths{ii} = [CFDfiblength DTIfiblength];
    FibPennation{ii} = [deepCFDpenn supCFDpenn deepDTIpenn supDTIpenn];
    
    prop30{ii} = p30; prop20{ii} = p20; prop10{ii} = p10;
    
    
end

cd('D:\Users\ghan078\Documents\corrected_data')

if ui == 2
    save('firstfour_FG.mat', 'CFDtracts','DTItracts','cfdOut','dtiOut','FibLengths','FibPennation','prop30','prop20','prop10');
elseif ui == 1
    save('firstfour_noFG.mat', 'CFDtracts','DTItracts','cfdOut','dtiOut','FibLengths','FibPennation','prop30','prop20','prop10');
end

% RMSE = sqrt(SSqE);
% 
% DataCells = cell(0);
% 
% % Filling in Row Labels
% Prefix = 'Subject0';
% for i = 2:5
%     DataCells{i,1} = [Prefix num2str(i-1)];
% end
% 
% % Filling in Column Labels
% Headers = [{'% within 10degrees'} {'% within 20degrees'} {'% within 30degrees'} ...
%     {'median deviation'} {'RMSE'} {'CFD fiber length'} {'DTI fiber length'} ...
%     {'fiber length RMSE'} {'CFD PCSA'} {'DTI PCSA'} ...
%     {'CFD pennation'} {'DTI pennation'} {'pennation RMSE'}];
% DataCells(1,2:14) = Headers;
% 
% %quick correction for assembling data
% % prop10(1:4) = []; prop20(1:4) = []; prop30(1:4) = [];
% % ThetaMed(1:4) = []; RMSE(1:4) = [];
% % CFDpcsa(1:4) = []; DTIpcsa(1:4) = []; medCFDpenn(1:4) = []; medDTIpenn(1:4) = []; RMSE_penn(1:4) = [];
% 
% % Assigning Data
% DataCells(2:5,2:6) = [prop10(1:4)' prop20(1:4)' prop30(1:4)' num2cell(ThetaMed(1:4))' num2cell(RMSE(1:4))'];
% DataCells(2:5,7:9) = [CFDmedflength(1:4)' DTImedflength(1:4)' RMSE_fl(1:4)'];
% DataCells(2:5,10:14) = [CFDpcsa(1:4)' DTIpcsa(1:4)' meddeepCFDpenn(1:4)' meddeepDTIpenn(1:4)' RMSE_penn_deep(1:4)'];
% 
% % Saving some stuff
% % cd ..
% % save('CFDDTIvars.mat', 'DataCells', 'CFDtracts', 'DTItracts', 'CFDlength',...
% %     'DTIlength','CFDpenn', 'DTIpenn', 'thetaOut', 'cfdOut', 'dtiOut');
% 
% %% Statistics
% 
% figure; hist(thetaOut{2,6}(:,1))
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w', 'facealpha',0.75);
% hold on
% hist(thetaOut{2,6}(:,2))
% h1 = findobj(gca,'Type','patch');
% set(h1,'facealpha',0.75);
% 
% 
% for i = 1:8
%     vecMat = [];
%     for j = 1:length(CFDtracts{1,i})
%         aVec = diff(CFDtracts{1,i}{j,1});
%         vecMat = [vecMat;mean(aVec,1)];
%     end
%     testvecs{i} = vecMat./repmat(sqrt(sum(vecMat.^2,2)),1,3);
%     %     testvecs{i} = unique(testvecs{i},'rows');
%     
%     dp = sum(testvecs{i} .* repmat([1,0,0],length(testvecs{i}),1),2);
%     ax_cfd{i} = unique(acosd( sign(dp).*dp ));
%     
% end
% 
% for i = 1:8
%     vecMat = [];
%     for j = 1:length(DTItracts{1,i})
%         aVec = diff(DTItracts{1,i}{j,1});
%         vecMat = [vecMat;mean(aVec,1)];
%     end
%     testvecs2{i} = vecMat./repmat(sqrt(sum(vecMat.^2,2)),1,3);
%     %     testvecs2{i} = unique(testvecs2{i},'rows');
%     
%     dp = sum(testvecs2{i} .* repmat([1,0,0],length(testvecs2{i}),1),2);
%     ax_dti{i} = unique(acosd( sign(dp).*dp ));
%     
% end
% 
% 
% figure; hist(ax_cfd{1,1})
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w', 'facealpha',0.75);
% hold on
% hist(ax_dti{1,1})
% h1 = findobj(gca,'Type','patch');
% set(h1,'facealpha',0.75);
% 
% figure; hist(CFDlength{1,4}(CFDlength{1,4}>20,:))
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w', 'facealpha',0.75);
% hold on
% hist(DTIlength{1,4})
% h1 = findobj(gca,'Type','patch');
% set(h1,'facealpha',0.75);
% 
% 
% x = [thetaOut{2,5}; thetaOut{2,6}; thetaOut{2,7}; thetaOut{2,8}];
% 
% y = [thetaOut{3,5}; thetaOut{3,6}; thetaOut{3,7}; thetaOut{3,8}];
% 
% z = [thetaOut{4,5}; thetaOut{4,6}; thetaOut{4,7}; thetaOut{4,8}];
% 
% [hx,px] = ttest(x(:,1),x(:,2));
% [hy,py] = ttest(y(:,1),y(:,2));
% [hz,pz] = ttest(z(:,1),z(:,2));
% 
% 
% 
% 
% 
% 
% 
% % x = [x ; thetaOut{2,6}; thetaOut{3,6}; thetaOut{4,6}];
% %
% % x = [x ; thetaOut{2,7}; thetaOut{3,7}; thetaOut{4,7}];
% %
% % x = real([x ; thetaOut{2,8}; thetaOut{3,8}; thetaOut{4,8}]);
% 
% 
%     
%     
%     %   % Old versions of file location
%     % CFDfolders = {'Male1_InvSeg_nostl_80percentProxApo\Design 1\Scenario 1\solver'...
%     %     'Female2_InvSeg_80percentProxApo\Design 1\Scenario 1\solver'...
%     %     'Male3_InvSeg_shaperedo_80percentProxApo_3\Design 1\Scenario 1\solver'...
%     %     'Female4_InvSeg_80percentProxApo\Design 1\Scenario 1\solver'};
% 
% 
%     
%             
%         
% %         idxnow = find(ismember(data.tracts_xyz',origins(j,:),'rows'));
% %         [a,b] = find(data.fibindex == idxnow);
% %         ai = data.fibindex(a,1);
% %         bi = data.fibindex(a,2);
% %         plot3(data.tracts_xyz(1,ai:bi),data.tracts_xyz(2,ai:bi),data.tracts_xyz(3,ai:bi))
% %         fibertracts_DTI{j} = [data.tracts_xyz(1,ai:bi)',data.tracts_xyz(2,ai:bi)',data.tracts_xyz(3,ai:bi)'];
% 
% 
