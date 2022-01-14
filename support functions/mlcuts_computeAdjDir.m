% This function is to compute adjacency matrix of graph based on direction 
% features in the following paper:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: direction: direction features in the current graph
%        sigma_T: similarity measure decay parameter on relative direction
%                 theta
% Output: adj_dir: adjacency matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adj_dir = mlcuts_computeAdjDir(direction,sigma_T)
% direction: (num_of_node,3);
%direction = nodeDir;
num_of_nodes = size(direction,1);
adj_dir= zeros(num_of_nodes);
for i = 1:num_of_nodes-1
    for j = i+1:num_of_nodes
        innervalue = inner(direction(i,:),direction(j,:));
        costheta = abs(innervalue)/(norm(direction(i,:),2)*norm(direction(j,:),2));
        
        % adjacent
        adj_dir(i,j) = exp(-(costheta-1).^2/sigma_T);% 0.04%costheta;% %costheta;%% 0.09 for sparse data%costheta; %%*costheta_mean2; 
        adj_dir(j,i) = adj_dir(i,j);
    end
end