% This supporting function for mLCuts to find the Gni for all the nodes in 
% the following paper. The thoery for 'point' or 'segment' are similar, but
% to clean the code, we separate the two functions.
%
% Wang J, Zhang M, Zhang J, Wang Y, Gahlmann A, Acton ST.
% Graph-theoretic post-processing of segmentation with application to dense
% biofilms. Submitted to IEEE TIP 2021.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: segLoc: super node location contains three kind of locaton:
%                segLoc{1,:}: all the voxels in all the segments;
%                segLoc{2,:}: centroids of the segments;
%                segLoc{3,:}: boundary points of the segments;    
%        hopNum: the number of hops in the neighborhood
% Output: Gni:  the multi-hop neighborhood for all the nodes
%               (the stored values are indices of nodes in the graph)
%         segmindist: pairwise seg shortest distance
%         adj_conn: adjacency matrix for pairwise connection between
%                   segments (1 for connected, 0 for no connection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Gni, segmindist,adj_conn] = findGn_segment(segLoc, hopNum,overlapArea)
%% Find Gn
num_of_nodes = size(segLoc,2);
% Compute connection adj matrix
adj_conn = zeros(num_of_nodes);
segmindist = zeros(num_of_nodes);
mindist_counts = zeros(num_of_nodes);
segDistLimit = sqrt(5);
% overlapArea = pi*4;

for i = 1:num_of_nodes-1
    for j = i+1:num_of_nodes
        [connect,mindist,mindistct] = Gn_checkSegConnection(i,j,segLoc,segDistLimit,overlapArea);
        if connect == 1
            adj_conn(i,j) = 1;
        end
        segmindist(i,j) = mindist;
        mindist_counts(i,j) = mindistct;
        segmindist(j,i) = segmindist(i,j);
        mindist_counts(j,i) =  mindist_counts(i,j);
        adj_conn(j,i) = adj_conn(i,j);
    end
end
% up to now, we have segdist, and adj_conn
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
