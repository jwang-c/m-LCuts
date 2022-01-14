% This function is the supporting function for post-segmentation using
% m-LCuts to rotate geometrical models(spherecylinders) to align with the
% identified cell directions after LCuts.
%
% For more information, please read/cite the following paper:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% [2] Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST and 
% Gahlmann A. Non-invasive single-cell morphometry in living bacterial 
% biofilms. Nature communications, 11(1), pp.1-13, 2020.
%
% Jie Wang, University of Virginia
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rot = postV_rotateBact(data) 
final_unit_vector = data(3:5);
rotation_vector = vrrotvec([0,0,1], final_unit_vector);
rotation_vector(4) = -rotation_vector(4);
Rot = vrrotvec2mat(rotation_vector);
end