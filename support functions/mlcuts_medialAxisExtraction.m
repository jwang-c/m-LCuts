% This supporting function for mLCuts to construct graph from 
% under-segmented masks.
%
% For more information, please refer to the following paper:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: Mask: the current segment
%        r_min: min radius for radius constrained medial axis extraction
%        r_max: max radius for radius constrained medial axis extraction
% Output: nodeLoc: coordinates of extracted point cloud (medial axes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodeLoc = mlcuts_medialAxisExtraction(Mask,r_min,r_max)
    seg2process = Mask;
    BW=bwdist(bwperim(seg2process));
    
    initialMask=zeros(size(seg2process));
    initialMask(seg2process==1)=1;  
    Mask = initialMask.*BW;

    %% limit radii range based on prior bacterial info;
    Skeleton = bwskel(Mask>r_min & Mask<=r_max);
    nodeLoc = [];
    for z = 1:size(Skeleton,3)
        [xstack,ystack] = find(Skeleton(:,:,z)==1);
        zstack = z*ones(size(xstack,1),1);
        nodeLoc = [nodeLoc;xstack,ystack,zstack];
    end
    %figure;plot3(nodeLoc(:,1),nodeLoc(:,2),nodeLoc(:,3),'.','MarkerSize',10);axis equal;