% This function is to find the direction of current component in the graph.
% Support function for finding principle direction in the following paper:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: currentSeg: current point cloud data
% Output: direction: principle direction vector,
%         eccentricity:  of current point cloud data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [direction,eccentricity] = Gn_pca(currentSeg)
Sub = ones(size(currentSeg));
Sub = mean(currentSeg,1).*Sub;
ShiftSeg = currentSeg - Sub;

CovM = cov(ShiftSeg);
[V,D] = eig(CovM);
DD = diag(D);
scale1 = DD(3);
if size(currentSeg,2) ==2
    x_v1= V(1,2);
    y_v1= V(2,2);
    direction = [x_v1,y_v1];
    eccentricity = sqrt(1-((1/sqrt(D(2,2))).^2/(1/sqrt(D(1,1))).^2));
else
    
    x_v1= V(1,3);
    y_v1= V(2,3);
    z_v1= V(3,3);
    direction = [x_v1,y_v1,z_v1];
    eccentricity = sqrt(1-((1/sqrt(D(3,3))).^2/(1/sqrt(D(2,2))).^2)); % 
end