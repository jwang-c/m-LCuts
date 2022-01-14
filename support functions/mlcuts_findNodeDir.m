% This supporting function for mLCuts to find all the node Dir features in:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: nodeLoc: node locations
%        mode: 1 - undersegment/point or 2- oversegment/segment
% Output: nodeLoc: coordinates of extracted point cloud (medial axes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodeDir, dist, adj_conn] = mlcuts_findNodeDir(nodeLoc, MODE)
%% Construct Gn in the paper
% Gni contain the neighborhood for all the nodes i
% ------------------------
% other parameters as mentioned in parameter validation section in the paper:
% these parameters can be estimated by different applications
hopNum = 5; % maximum hop-level number % >5
overlapArea = 3*pi; % minimum neighboring segment contacting area (for segment mode)
% ------------------------
if MODE == 1
    [Gni, dist, adj_conn] = findGn_point(nodeLoc, hopNum); 
    nodeDir = zeros(size(nodeLoc,1),3);
elseif MODE == 2
    [Gni, dist,adj_conn] = findGn_segment(nodeLoc, hopNum,overlapArea);
    nodeDir = zeros(size(nodeLoc,2),3);
end
% Note: to avoid duplicated computation of nodeLoc distances and connections
% we compute dist and adj_conn ahead of time and pass them to the later
% functions. They are not shown in the algorithm shown in the paper.
%% Find major direction
%nodeDir = zeros(size(nodeLoc,1),3);
%nodeLocs = [];
%for a = 1:size(nodeLoc,2)
%    nodeLocs = [nodeLocs;nodeLoc{2,a}];
%end
%figure;axis equal;
for i = 1:size(nodeDir,1)
    nodeDir(i,:) = mlcuts_findMajorDir(i,Gni,nodeLoc,MODE, hopNum, dist);
    %hold on;
    %plot3([0+nodeLocs(i,1), nodeDir(i,1)*4+nodeLocs(i,1)],[0+nodeLocs(i,2),nodeDir(i,2)*4+nodeLocs(i,2)],[0+nodeLocs(i,3), nodeDir(i,3)*4+nodeLocs(i,3)],'r');
end