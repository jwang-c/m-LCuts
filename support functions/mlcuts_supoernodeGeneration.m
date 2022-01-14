% This supporting function for mLCuts in graph construction in the
% supernode/segment/oversegmentation mode.
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: Mask: current oversegmented regions
% Output: nodeLoc: all coordinates of supernodes in column 1
%                  centroids locations of supernodes in column 2
%                  boundary points of supernodes in column 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nodeLoc = mlcuts_supoernodeGeneration(Mask)
    %% node locations
    clear nodeLoc; 
    %figure;isosurface(Mask,0)
    for n = 1:max(Mask(:))
        coord = [];
        for z = 1:size(Mask,3)
            currentseg = Mask(:,:,z)==n;
            if sum(currentseg(:))>0            
                [row,col] = find(currentseg);
                stack = z*ones(size(row,1),1);
                coord = [coord;row col stack];
            end
        end
        nodeLoc{1,n} = coord; 
        nodeLoc{2,n}= mean(coord,1);
        %figure;plot3(coord(:,1),coord(:,2),coord(:,3),'.');hold on;axis equal;axis off;
    end

    %% find boundary points to calculate distance
    for i = 1:size(nodeLoc,2)
        x = nodeLoc{1,i}(:,1);
        y = nodeLoc{1,i}(:,2);
        z = nodeLoc{1,i}(:,3);
        currentComp = zeros(size(Mask));
        for xx = 1: length(x)
            currentComp(x(xx),y(xx),z(xx)) = 1;
        end
        zBoundaries = [];
        for z = 1:size(Mask,3)
            B= bwboundaries(currentComp(:,:,z));
            if size(B,1)~=0
               zBoundaries = [zBoundaries;B{1,1} z*ones(size(B{1,1},1),1)];
            end
        end
        for x = 1:size(Mask,1)
            B= bwboundaries(squeeze(currentComp(x,:,:)));
            if size(B,1)~=0
               zBoundaries = [zBoundaries;x*ones(size(B{1,1},1),1) B{1,1}];
            end
        end
     %         zBoundaries
        nodeLoc{3,i} = unique(zBoundaries,'rows');
    end