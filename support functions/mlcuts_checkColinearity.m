% This function is the stopping criterion step in the following paper:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: sinComp: current single component
%        nodeLoc: node/supernode location
%        adj_conn: adjacency matrix for connections in graph
%        MODE: current running mode
% Output: sinCompLinear: '1' if the current compoenent satisfied stopping criterion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sinCompLinear = mlcuts_checkColinearity(sinComp,nodeLoc,adj_conn,MODE)

%compute the current principle direction and eccentricity of the current
%component
if MODE == 1
    currentNodes = nodeLoc(sinComp,:);  
elseif MODE == 2
    currentLocs = [];
    for s = 1:size(sinComp,1)
        currentLocs = [currentLocs;nodeLoc{1,sinComp(s)}];% size of current seg
    end 
    % represent each segment roughly by the centroids
    nodeCentroids = [];
    for c = 1:size(nodeLoc,2)
        nodeCentroids = [nodeCentroids;nodeLoc{2,c}];
    end
    currentNodes = nodeCentroids(sinComp,:);
end
[sinCompDirC,ecc] = Gn_pca(currentNodes);

%% check variations
if size(currentNodes,2) == 2 % 2D data case
    translate = mean(currentNodes,1);
    normal = [-sinCompDirC(2),sinCompDirC(1)];
    d = zeros(size(currentNodes,1),1);
    for c = 1: size(currentNodes,1)
        vector = currentNodes(c,:)-translate;
        d(c) = inner(normal,vector) / norm(normal,2);
    end
elseif size(currentNodes,2) == 3
    Sub = ones(size(currentNodes));
    translate = mean(currentNodes,1).*Sub; 
    currentPoints = currentNodes-translate;
    line = sinCompDirC;
    d = zeros(size(currentNodes,1),1);
    for c = 1: size(currentPoints,1)
        x0= currentPoints(c,:);
        d(c) = norm(cross(x0,line),2) / norm(line,2);
    end
end
%% avoid big gap and they are in the same group
% Both methods work, now using the adj_conn to compute
% ordermin = sort(dist(sinComp,sinComp));
% checkmin = ordermin(2,:);
% bigGap = sum(checkmin>=sqrt(5))>0;
checkconn = sum(adj_conn(sinComp,sinComp),1);
bigGap = sum(checkconn<=1)>2;
%% use convexhull
% mode =1 means seg, 0 for points
if MODE ==2 % segment
    [k,vol] = convhulln(currentLocs);
    s = vol/size(currentLocs,1);
else
    s = 0;
end

%% stopping criterion
if sum(d>4)<=1  && bigGap == 0 && s < 1.2 % s<1.8
    sinCompLinear = 1;
else 
    sinCompLinear = 0;
end
    