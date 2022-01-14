% This main function for mLCuts for data clustering with the application to 
% post-processing dense biofilms.
%
% For more information, please refer to the following papers:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Other reference: 
% [2] Wang J, Batabyal T, Zhang M, Zhang J, Aziz A, Gahlmann A, Acton ST.
% Lcuts: Linear Clustering of Bacteria Using Recursive Graph Cuts. In 2019
% IEEE International Conference on Image Processing (ICIP) 2019 Sep 22 (pp.
% 1575-1579). IEEE.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: Mask: region of interest;
%        MODE: "under" or "over".
% Output: final_Groups: The final clustered groups;
%         postV:  (Optional) reconstructed volume
% Parameters: sizeLimit,distLimit,σD,σT
% Other inputs for reconstruction and result save: 
%       expname (experiment name), L (initial segmentation result),
%       need_splitn (the current index for ROI), showI, saveI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [final_Groups,postV] = mLCuts(Mask,MODE,sizeLimit,distLimit,r,sigmaD,sigmaT,expname,L,need_splitn,showI,saveI)%, GT,tG,roi)

%% Find nodes with features:
if MODE == 1
    nodeLoc = mlcuts_medialAxisExtraction(Mask,r(1),r(2));
    nodeLoc = mlcuts_cleanNode(nodeLoc); %(selective) just in case there are nodes that do not have neighbors -- clear them
    [nodeDir, dist, adj_conn] = mlcuts_findNodeDir(nodeLoc, MODE); % find all nodeDir with algorithm 2 in the paper
elseif MODE == 2
    nodeLoc = mlcuts_supoernodeGeneration(Mask);
    [nodeDir, dist, adj_conn] = mlcuts_findNodeDir(nodeLoc, MODE); % note the dist for supernode is the min dist between two segment
end
%% Construct the graph using eq.(1) in the paper [1]
WD = exp(-(dist).^2./(sigmaD^2)); %eq.(2)
WT = mlcuts_computeAdjDir(nodeDir, sigmaT); %eq.(5)
W = WD.*WT; % eq.(1), the penalty(eq.(7)) will be calculated inside of each Cut
%% Recursive cut: algorithm 3 in paper [1]
% -------initialize----------
clear groups;
if MODE == 1
    num_of_nodes = size(nodeLoc,1);
elseif MODE== 2
    num_of_nodes = size(nodeLoc,2);
end
nodes = 1:1:num_of_nodes;
N = nodes.';
it = 1;  % this is n in the paper [1]
status = 0;
groups{it} = [];
sinComp = N;
adj = W;
r = (r(2)+r(1))/2;
% -------compute advanced LCuts------
[groups,it] = mlcuts_LCuts(groups,it,sinComp,status,nodeLoc,nodeDir,adj,adj_conn,sizeLimit,distLimit,r,MODE);

%% show results and reconstruct the biofilm
clear final_Groups;
% Note: groups are indices, final_groups are point locations
if MODE == 1

    V_cutted = zeros(size(L)); % to save in the image
    count = 1;
    for i = 1:size(groups,2)  
        curComps = [];
        if size(groups{1,i},1) > 0
            for j = 1:size(groups{1,i},1)
                curComps = [curComps; nodeLoc(groups{1,i}(j),:)];
                V_cutted(nodeLoc(groups{1,i}(j),1),nodeLoc(groups{1,i}(j),2),nodeLoc(groups{1,i}(j),3))=count;
            end
            final_Groups{count} =curComps;
            count = count+1;
        end
    end

elseif MODE == 2
    for i = 1:size(groups,2)

        curComps = [];
        for j = 1:size(groups{1,i},1)
            curComps = [curComps; nodeLoc{1,groups{1,i}(j)}];
        end
        final_Groups{i} =curComps;
    end 

end
%% show final groups
colorm = [];
for i = 1:size(final_Groups,2)
    colorm = [colorm;rand(1,3)];
end

if showI
    figure(10);
    hold on;axis equal;
    cmap = jet(size(final_Groups,2));
    for i = 1:size(final_Groups,2)
        cm = colorm(i,:);
        y = final_Groups{1,i}(:,1);
        x = final_Groups{1,i}(:,2);
        z = final_Groups{1,i}(:,3);
        plot3(x,y,z,'.','Color',cm ,'MarkerSize',20);axis equal off;view([-37.5,30])
        hold on;
    end
end

if saveI == 1
    filename = [expname+'_mlcutsresult.png'];
    saveas(gcf,filename);
    filename = [expname+'_mlcutsresult.fig'];
    saveas(gcf,filename);
end

%% model reconstruction and genarate post ROI volume
%Reconstruct the biofilm by fitting geometrical models
if MODE == 1
    seg2process = L == need_splitn;
    %figure;isosurface(seg2process,0);axis equal;
    clear final_Comp post_radii;
    for count =  1:size(final_Groups,2)
        [BW_medial,~]=bwdist(bwperim(seg2process));
        post_radii{1,count} = postV_findMaxRadii(final_Groups{1,count},BW_medial);

    end
    final_Comp = postV_LCuts2Surfaces(final_Groups,post_radii);
    label = 1;
    postV = zeros(size(L));
    volume_lowerbound = 4/3*pi*4*4*4;
    for b = 1:size(final_Comp,2)
        current_Bact = round(final_Comp{1,b});
        current_Bact = unique(current_Bact,'row');
        current_Bact(current_Bact<=0)=1;
        current_Bact(current_Bact(:,3)>=size(L,3),3)=size(L,3);
        current_Bact(current_Bact(:,2)>=size(L,2),2)=size(L,2);
        current_Bact(current_Bact(:,1)>=size(L,1),1)=size(L,1);
        if size(current_Bact,1)> volume_lowerbound
            for x = 1:size(current_Bact,1)
                postV(current_Bact(x,1),current_Bact(x,2),current_Bact(x,3)) = label;
            end
            label = label+1;
        end
    end
    %%%%%%%%%
    if showI
        figure;
        hold on;axis equal;
        for i = 1:size(final_Comp,2)
            cm = colorm(i,:);
            y = final_Comp{1,i}(:,1);
            x = final_Comp{1,i}(:,2);
            z = final_Comp{1,i}(:,3);
            plot3(x,y,z,'.','Color',cm ,'MarkerSize',10);axis equal off;view([-37.5,30])
            hold on;
        end
    end
    %
    if saveI == 1
        filename = [expname+'_mlcuts_reconstruction.png'];
        saveas(gcf,filename);
        filename = [expname+'_mlcuts_reconstruction.fig'];
        saveas(gcf,filename);
    end
elseif MODE == 2
    postV = zeros(size(L));
    for b = 1:size(final_Groups,2)
        current_Bact = final_Groups{1,b};
        for x = 1:size(current_Bact,1)
            postV(current_Bact(x,1),current_Bact(x,2),current_Bact(x,3)) = b;
        end 
    end
end

%% save the workspace
%savename = expname+'.mat';
%save(savename);