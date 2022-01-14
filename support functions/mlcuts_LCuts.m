% This is the main function to run recursive LCuts -- algorithm 3 in the
% following paper:
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
%        The current number of groups it in the data;
%        The status of recursion;
%        The current groups in the data.
% Output: The updated groups in the data;
%         The updated number it of groups.
% Parameters: sizeLimit,distLimit,σD,σT,MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
function [groups,it] = mlcuts_LCuts(groups,it,sinComp,status,nodeLoc,nodeDir,adj,adj_conn,sizeLimit,distLimit,r,MODE)
    % ----- update the graph ------ %
    if MODE == 1 % points
        currentSize = size(sinComp,1); 
    elseif MODE == 2 % segment
        currentSize = 0;
        for s = 1:size(sinComp,1)
            currentSize = currentSize+size(nodeLoc{1,sinComp(s)},1);% size of current seg, for seg nodeLoc is a 1*3 cell
        end 
    end
    % ----- check recursion stopping criterion ------ %
    if currentSize <= sizeLimit(1) && currentSize > sizeLimit(2) && size(sinComp,1)>=2
        sinCompLinear =  mlcuts_checkColinearity(sinComp,nodeLoc,adj_conn,MODE);
    else
        sinCompLinear = 0;
    end  

    if status == 1 || sinCompLinear ==1 || currentSize <= sizeLimit(2) || size(sinComp,1)< 2
        groups{it} = sinComp;
        it = it+1;            
    else
       % ----- bi-partition ------ %
       Comps = mlcuts_Cut(sinComp,nodeLoc,nodeDir,adj,r,distLimit,MODE);
       if size(Comps,2)==1 || isempty(Comps{1}) || isempty(Comps{2}) % complete
           status = 1;
            [groups,it] = mlcuts_LCuts(groups,it,sinComp,status,nodeLoc,nodeDir,adj,adj_conn,sizeLimit,distLimit,r,MODE);
       else % continue
           status = 0;
           sinComp = Comps{1};
           [groups,it] = mlcuts_LCuts(groups,it,sinComp,status,nodeLoc,nodeDir,adj,adj_conn,sizeLimit,distLimit,r,MODE);
           sinComp = Comps{2};
           [groups,it] = mlcuts_LCuts(groups,it,sinComp,status,nodeLoc,nodeDir,adj,adj_conn,sizeLimit,distLimit,r,MODE);
       end           
    end

   return
end