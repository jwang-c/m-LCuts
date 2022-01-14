% This supporting function for mLCuts to find the Gni for all the nodes in 
% the following paper. The thoery for 'point' or 'segment' are similar, but
% to clean the code, we separate the two functions.
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: nodeLoc: node location
%        hopNum: the number of hops in the neighborhood
% Output: Gni:  the multi-hop neighborhood for all the nodes
%               (the stored values are indices of nodes in the graph)
%         dist: pairwise distance between nodes
%         adj_conn: adjacency matrix for pairwise connection between nodes
%         (1 for connected, 0 for no connection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Gni, dist, adj_conn] = findGn_point(nodeLoc, hopNum)
%% Find Gn
num_of_nodes = size(nodeLoc,1);
% Compute connection adj matrix
adj_conn = zeros(num_of_nodes);
dist = zeros(num_of_nodes);

dT = 2; % Thresholding value to consider a neighborhood in the paper

for i = 1:num_of_nodes-1
    for j = i+1:num_of_nodes
        dist(i,j) = norm(nodeLoc(i,:)-nodeLoc(j,:),2);
        dist(j,i) = dist(i,j);
        if dist(i,j) < dT
            adj_conn(i,j) = 1;
            adj_conn(j,i) = adj_conn(i,j);
        end
    end
end
% up to now, we have dist, and adj_conn
% other stats:
% temp_dist = dist;
% temp_dist(dist == 0)=max(dist(:));
% mindist = min(temp_dist,[],2);
% vardist = sum((mindist-mean(mindist)).^2)/length(mindist);
%% Define Gni for all the nodes:
% note: here Gni is the 5-hop neighborhood coonnection list
clear Gni; 
%hopNum = 5;
for n = 1:num_of_nodes
    hop = 1;
    currentConn = n;
    preConn = [];
    while hop <= hopNum
        checkConn = currentConn;
        currentConn = [];
        for num = 1:length(checkConn)
            Conn = Gn_checkConnection(checkConn(num),adj_conn,preConn);
            currentConn = [currentConn,Conn];
            currentConn = unique(currentConn);
        end
        Gni{n,hop} = currentConn; % Gni
        preConn = [preConn,currentConn];
        hop = hop+1;
    end
end
