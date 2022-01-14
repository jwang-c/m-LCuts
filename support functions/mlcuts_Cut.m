% This is the function to compute bi-partition solution using N-Cut in the
% following papers:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% Other references:
%
% [2] Wang J, Batabyal T, Zhang M, Zhang J, Aziz A, Gahlmann A, Acton ST.
% Lcuts: Linear Clustering of Bacteria Using Recursive Graph Cuts. In 2019
% IEEE International Conference on Image Processing (ICIP) 2019 Sep 22 (pp.
% 1575-1579). IEEE.
%
% Jie Wang, University of Virginia, VIVA lab
% Last update: May-28-2021 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: The current nodes sinComp in the single component;
%        The node features nodeLoc and nodeDir; 
%        The current adjacency matrix adj of sinComp;
%        The estimated radius of cells, r;
%        The distance limit between two nodes in the same group;
% Output: The updated groups in the data;
%         The updated number it of groups.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Comps = mlcuts_Cut(sinComp,nodeLoc,nodeDir,adj,r,distLimit,mode)
nodes = sinComp;
num_of_nodes = size(sinComp,1);
%% update current adj without penalty
currentAdj = adj(nodes,nodes);

%% compute penalty, equation (7) in paper [1]:
if mode == 2 % use centroids to indicate each supernode
    nodeCentroids = [];
    for s = 1:size(nodeLoc,2)
        nodeCentroids = [nodeCentroids;nodeLoc{2,s}];
    end
    Locs = nodeCentroids(nodes,:);
    %Locs = segLoc{3,nodes};
    delta = zeros(num_of_nodes);
    for i = 1:num_of_nodes
        for j = i+1:num_of_nodes       
            dist = norm(Locs(i,:)-Locs(j,:),2);
            if dist <= distLimit
                % compare dij and dji select the large one
                % point j to dir i
                translatei = Locs(i,:); 
                currentLocj = Locs(j,:)-translatei;
                linei = nodeDir(nodes(i),:);
                dji= norm(cross(currentLocj,linei),2) / norm(linei,2);
                % point j to dir i
                translatej = Locs(j,:); 
                currentLoci = Locs(i,:)-translatej;
                linej = nodeDir(nodes(j),:);
                dij= norm(cross(currentLoci,linej),2) / norm(linej,2);

                distij = max(dij,dji);
                delta(i,j) = exp((-distij^2)/((r+1)^2));
            else
                delta(i,j) = 0;
            end
            delta(j,i) = delta(i,j);
        end
    end
else
    delta = 1;
end
   
%% update the whole adjacency matrix
adj = currentAdj.*delta;%%;
deg = diag(sum(adj,1)); % D in the paper
invdeg = sqrt(diag(1./sum(adj,1))); % D^-0.5
symL = invdeg*(deg-adj)*invdeg;% symmetric laplacian
% check if there are rows all equals to nan/inf
if any(sum(symL==Inf,2)== num_of_nodes)
    status = zeros(size(symL,1),1);
    status((sum(symL==Inf,2)== num_of_nodes))=1;
elseif any(sum(isnan(symL),2)== num_of_nodes)
    status = zeros(size(symL,1),1);
    status((sum(isnan(symL),2)== num_of_nodes))=1;
else %nomralized cut
    [UU,VV] = eig(symL);

    sVV = sort(diag(VV));
    indF = find(diag(VV)==sVV(2));
    status = zeros(size(adj,1),1);

    for k = 1:size(adj,1)
       if UU(k,indF) < 0
           status(k) =0;
       else
           status(k)=1;
       end
    end
end
% 
ind1 = find(status);
Comps{1} = nodes(ind1,: );
ind0 = find(1-status);
Comps{2} = nodes(ind0,: );
