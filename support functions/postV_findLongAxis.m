% This is the support function for reconstruct biofilms --  
% reconstruct biofilm by fitting spherocylinder model
%
% For more information, please refer to the following papers:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Other reference: 
% [2] Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST and 
% Gahlmann A. Non-invasive single-cell morphometry in living bacterial 
% biofilms. Nature communications, 11(1), pp.1-13, 2020.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function longaxis = postV_findLongAxis(currentBact)
    dist = zeros(size(currentBact,1));
    for c = 1:size(currentBact,1)
        for r = 1:size(currentBact,1)
            dist(c,r) = norm(currentBact(c,:)-currentBact(r,:),2);
        end
    end
    
    longaxis = max(dist(:));