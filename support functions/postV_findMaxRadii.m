% This is the support function for reconstruct biofilms in mLCuts --  
% find maximum radii along medial axes.
%
% For more information, please refer to the following papers:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% [2] Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST and 
% Gahlmann A. Non-invasive single-cell morphometry in living bacterial 
% biofilms. Nature communications, 11(1), pp.1-13, 2020.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radii = postV_findMaxRadii(current_postseg,BW_medial)

L = size(current_postseg,1);
radii = [];
for dl = 1:L
    currentradii = BW_medial(current_postseg(dl,1),current_postseg(dl,2),current_postseg(dl,3));
    radii = [radii; currentradii];
end
    
