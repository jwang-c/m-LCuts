% This is the support function for reconstruct biofilms in mLCuts --  
% reconstruct biofilm by fitting spherocylinder model
%
% For more information, please refer to the following papers:
%
% [1] Wang J, Zhang M, Zhang J, Wang Y, Gahlmann A, Acton ST.
% Graph-theoretic post-processing of segmentation with application to dense
% biofilms. Submitted to IEEE TIP 2021.
%
% [2] Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST and 
% Gahlmann A. Non-invasive single-cell morphometry in living bacterial 
% biofilms. Nature communications, 11(1), pp.1-13, 2020.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [currentBact,Rot,centerMove] = postV_createModel(data,r)
    % if component length is too short, use half the radius
    if data(1) < 8.5
        r = (0.5*(data(2)));
    end
    l = (data(1)-2*r); % the length of the cylinder part
    if l < 1
        l=1;
    end
    grid_size = 0.5;
    ellipsoid_template = postV_spherocylinder_uniformPoints(l,r,grid_size);
    Rot = postV_rotateBact(data);
    currentCentroid = data(6:8);
    centerMove = (currentCentroid).*ones(size(ellipsoid_template,1),3);
	currentBact = ellipsoid_template * Rot + centerMove;
end
