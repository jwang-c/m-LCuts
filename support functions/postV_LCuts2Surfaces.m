% This is the support function for reconstruct biofilms in mLCuts --  
% reconstruct biofilm by fitting spherocylinder model
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
function Bact = postV_LCuts2Surfaces(post_segments,post_radii)

%% for each segment after LCuts, find direction and length
stats = zeros(size(post_segments,2),8);
for num = 1:size(post_segments,2)
    currentBact = post_segments{1,num};   
    if size(currentBact,1)>1 % in case there is only one point in the group
        r_min = min(post_radii{1,num});
        stats(num,1) = postV_findLongAxis(currentBact)+2*r_min; %length, note that each segment is shorter than the actual cell length by r_min*2 due to medial axis extraction, so we need to add this back
        stats(num,2) = max(post_radii{1,num})+1; %radius
        [direction,~] = Gn_pca(currentBact);
        stats(num,3:5) = direction; %direction
        stats(num,6:8) = mean(currentBact,1); %center
    else
        stats(num,2) = -1;
    end
end
statistics = stats(stats(:,1)<=90 & stats(:,2)>0,:);

%% fit with spherecylinders based on the statistics
for i = 1:size(statistics,1)
    r = statistics(i,2);
    [Bact{i},Rot{i},centerMove{i}] = postV_createModel(statistics(i,:),r);
%     figure(2);hold on;
%     plot3(round(Bact{1,i}(:,2)),round(Bact{1,i}(:,1)), round(Bact{1,i}(:,3)),'.');  axis equal;title('Reconstructed cells by fitting models after LCuts');
end



