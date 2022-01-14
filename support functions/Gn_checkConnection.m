% This supporting function for graph construction in mLCuts to check what 
% are the neighboring nodes in the neighborhood for the current node n.
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: n: current node n
%        adj_conn: node connection adjacency matrix
%        preConn: the previous neighbors already found in the neighorhod Gn
% Output: conn: new neighbors (indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conn = Gn_checkConnection(n,adj_conn,preConn)

nidx = find(adj_conn(n,:));

Lia = ismember(nidx,preConn);
conn = nidx(~Lia);