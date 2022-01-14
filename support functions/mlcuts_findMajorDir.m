% This supporting function for mLCuts to find the major direction feature
% using Algorithm 2 in the following paper:
%
% [1] Wang J, Zhang M, Zhang J, Wang Y, Gahlmann A, Acton ST.
% Graph-theoretic post-processing of segmentation with application to dense
% biofilms. Submitted to IEEE TIP 2021.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: i: the index for the current node
%        Gni: multi-hop neighborhood connection list
%        nodeLoc: all the locations of nodes
%        mode: '1':undersegment/point or '2':oversegment/segment
%        ---- some extra input pass from the former algorithms to avoid duplicated computation ----
%        hopNum: the number of hops in Gni (the largest hop number)
%        dist: each attribute is a shortest distance between two nodes. In
%              the case of super node, it is the shortest distance between
%              all the points in the two segment.
%        adj_conn: adjacency matrix for connection in the graph Gn
%        -------------------------------------------------------------------------------------------
% Output: nodeDir: node major direction feature for the current node i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodeDir = mlcuts_findMajorDir(i,Gni,nodeLoc, MODE, hopNum, dist)
%% ---------------------- mode for point/undersegment ---------------------
% line 2 -- 5 of Algorithm in paper [1]
if MODE == 1
    %% Preparation: 
    % find the current Gni (multi-hop neighborhood), clear the 
    % duplicated points
    currentPoint = nodeLoc(i,:);
    connNodes = []; % indices of current neighboring node to currentPoint i
    for c = 1:hopNum
        connNodes = [connNodes Gni{i,c}];
    end
    selectIdx = unique(connNodes); % remove the duplicated neighbors
    for s = 1: length(selectIdx)
       if selectIdx(s)~= i
          currentPoint = [currentPoint; nodeLoc(selectIdx(s),:)];
       end
    end
    %--------------------------------------------------------------------------
    %% Find candidate directions: 
    % sort the distance between current node (currentPoint(1,:)) and other nodes
    for p = 1 : size(currentPoint,1)
        distPoint(p) = norm(currentPoint(p,:)-currentPoint(1,:),2);
    end
    [~,rankedIdx] = sort(distPoint);
    currentNode = currentPoint(rankedIdx,:);

    % find the candidate orientations in the Gni
    candidate = zeros(size(currentNode,1)-1+size(currentNode,1)-2,size(currentNode,2));
    for j = 1:size(currentNode,1)-1
        candidate(j,:) = currentNode(j+1,:)-currentNode(1,:);
    end
    % Adding more pairwise directions from the closest node to current node i can
    % help to provide more stable direction (selective)
    for jj = 1:size(currentNode,1)-2
        candidate(size(currentNode,1)-1+jj,:) = currentNode(jj+2,:)-currentNode(2,:);
    end
%% ----------------------- mode for segment/oversegment -------------------
% line 6 -- 9 of Algorithm in paper [1]
elseif MODE == 2
    %% Preparation: 
    % find the current Gni (multi-hop neighborhood), clear the 
    % duplicated points
    connNodes = i;
    for c = 1:size(Gni,2)
        connNodes = [connNodes Gni{i,c}];
    end
    selectIdx = unique(connNodes);

    segLoc = nodeLoc; % in the case of segment, the nodeLoc include three different styles
    %% Find candidate directions:
    if size(selectIdx,2) >1
        if size(segLoc{1,i},1) >=  4/3*pi*4*4*8% single compoment is long enough to find the dominate direction
            candidate = Gn_pca(segLoc{1,i});
        else
            %% Find candidate directions:
            % sort the distance between current node (currentSeg(1,:)) and other nodes
            dists = dist(i,selectIdx(1:length(selectIdx)));
            [~,rankedIdx] = sort(dists);
            currentNode = selectIdx(rankedIdx); % selected supernode neighbors in Gni
            %% find the orientations 
            candidate = zeros(size(currentNode,2)-1,3);%+size(sPoints,2)-2
            for j = 1:size(currentNode,2)-1
                currentSeg = [segLoc{1,currentNode(1)};segLoc{1,currentNode(j+1)}];
                candidate(j,:) = Gn_pca(currentSeg);
            end
        end
    else % if there's no neighbors other than the current segment
       candidate = Gn_pca(segLoc{1,selectIdx}); 
    end

end
%--------------------------------------------------------------------------
if size(candidate,1) == 1 % only one candidate
    nodeDir = candidate;
else
    %% line 12-16 in algorithm 2 in the paper [1]
    % compute pairwise direction similarity between two directions in
    % candidate using eq.(4)
    S = zeros(size(candidate,1)); 
    for o = 1: size(candidate,1)
        for oo = o+1: size(candidate,1)
            innervalue = inner(candidate(o,:),candidate(oo,:));
            S(o,oo) = abs(innervalue)/(norm(candidate(o,:),2)*norm(candidate(oo,:),2));
            S(oo,o) = S(o,oo);
        end
    end
    Smax = max(S,[],1);
    removeIdx = find(Smax<median(Smax));
    outlier = removeIdx(removeIdx<=size(currentNode,1)-1)+1;
    %outlier = removeIdx+1;
    %remove4sure = removeIdx(removeIdx<=size(candidate,1)-1)+1;
    clearedPoints = currentNode;
    clearedPoints(outlier,:)=[];
    Gni_minus=[];  % Gni- in the paper: outlier-removed neighborhood
    if MODE == 1
        Gni_minus = clearedPoints;
    elseif MODE == 2      
        for ns = 1:size(clearedPoints,2)
            Gni_minus = [Gni_minus;segLoc{1,clearedPoints(ns)}];
        end
    end
    if size(Gni_minus,1)<=1 % in case there is only one point in the graph
        Gni_minus = currentNode;
    end
    nodeDir = Gn_pca(Gni_minus); 
end
