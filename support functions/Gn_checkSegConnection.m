% This function will check if segment i is the neighbor of segment j by
% checking the shortest distance between two segment. In the code below, we
% use the boundary points on the segment to make the process faster.
%
% Support function for finding Gn in the following paper:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: i,j: the indices for two segment
%        segLoc: supernode/segment location
%        dT: shortest distance to consider a neighbor/connection
%        preConn: the previous neighbors already found in the neighorhod Gn
% Output: conn: new neighbors (indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [connect,mindist,mindist_count] = Gn_checkSegConnection(i,j,segLoc,dT,overlapArea)
    size_i = size(segLoc{1,i},1);
    size_j = size(segLoc{1,j},1);
    count = 0;
    dists = 100*ones(size_i,size_j); % initialize with a large number so that segments that are apart farther than 
    for a = 1:size_i
        for b = 1:size_j
            dist = norm(segLoc{1,i}(a,:)-segLoc{1,j}(b,:),2);
            if dist < dT
                count = count + 1;
            end
            dists(a,b) = dist;
        end
    end
    if count >= overlapArea%4
        connect = 1;
    else
        connect = 0;
    end
    mindist = min(dists(:));
    mindist_count = count;
%end