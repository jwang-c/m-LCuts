% This supporting function for mLCuts to clean the noisy nodes if they do 
% not have any neighbors.
% This part is mainly used for the original LCuts [2], where the point cloud 
% data is more noisy.
%
% Reference:
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% [2] Wang J, Batabyal T, Zhang M, Zhang J, Aziz A, Gahlmann A, Acton ST.
% Lcuts: Linear Clustering of Bacteria Using Recursive Graph Cuts. In 2019
% IEEE International Conference on Image Processing (ICIP) 2019 Sep 22 (pp.
% 1575-1579). IEEE.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: inputdata_o: original nodeLoc
% Output: nodeLoc: cleaned nodeLoc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodeLoc = mlcuts_cleanNode(inputdata_o)
%% clear the nodes: delete the data that has no neighbors within a distance

uinput = double(unique(inputdata_o,'rows'));
pointDistLimit = 5; % this is based on how to determine if there are noisy data
num_of_nodes = size(uinput,1);
checkdist = zeros(num_of_nodes);
for i = 1:num_of_nodes-1
    for j = i+1:num_of_nodes
        checkdist(i,j) = norm(uinput(i,:)-uinput(j,:),2);
        checkdist(j,i) = checkdist(i,j);
    end
end
check = checkdist<pointDistLimit;
checkNodes = sum(check,2);
inputdata = uinput(checkNodes~=1,:);
nodeLoc = inputdata;
%figure;plot3(nodeLoc(:,1),nodeLoc(:,2),nodeLoc(:,3),'.');axis equal;axis off;
