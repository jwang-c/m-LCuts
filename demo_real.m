% This is the demo code for experiments on real data follow the two-step 
% strategy in the paper:
%
% [1] J. Wang, M. Zhang, J. Zhang, Y. Wang, Andreas Gahlmann, and S. T. Acton,
% “Graph-theoretic Post-processing of Segmentation with Application to Dense 
% Biofilms.” IEEE Transaction on Image Processing, 30, 8580-8594.(2021)
%
% The post-processing needs may be uncatogorized. See section IV-E in [1].
% First run mode 1: undersegmentation mode; 
% Then run mode 2: oversegmentation mode.
%
% Other related papers:
% [2] BCM3D, DOI: 10.5281/zenodo.4088658
% [3] LCuts, DOI: 10.1109/ICIP.2019.8803064
%
% Latest update: Jan 14, 2021, Jie Wang
%%
clc;clear all; close all;
addpath('support functions');
%% Load data that requires post-processing
% The current demo image is tiff file.
load data_frame10_seg.mat
%%%%% Demo image is labeled segmentation results from [2] %%%%%

%% Sequential step 1: find undersegmented clusters that needs further splitting
% The current selection includes the clusters that is larger than a regular
% bacterial size, or the ones do not present a collinear cell shape.

%%% Parameters %%% 
% based on image resolution and bacterial size
VOLUME_UPPERBOUND = 4/3*pi*4*16+pi*16*50;
CONCAVITY_RATIO = 1.8;
%%%%%%%%%%%%%%%%%%

need_split = [];
labels = sort(unique(L(:)));
num_e = size(unique(L(:)),1)-1;
for n = 1:num_e
    temp_showV = zeros(size(L)); % an intermediate parameter to plot the volume
    temp_showV(L==labels(n+1)) = n;
    figure(1);isosurface(temp_showV,n-1);colormap colorcube;axis equal;axis off;title('Original input volume');
    cc = regionprops3((temp_showV==n),'ConvexVolume');
    s = cc.ConvexVolume/sum(temp_showV==n,'all');
    if size(s,1) >1 
        s = 0;
    end 
    if sum(temp_showV==n,'all')> VOLUME_UPPERBOUND || s>CONCAVITY_RATIO % size and concavity check
       need_split = [need_split;labels(n+1)]; % idx for clusters needs post-processing
    end
end

%% Run m-LCuts mode 1
close all;
%%% Parameters and initializations %%%
showI = 0;
saveI = 0;
distLimit = 5;
sigmaD = 2;
sigmaT = 0.1;
r = [2,6];
sizeLimit = [40,8];
MODE = 1;
postV1 = zeros(size(L));  % ROI volumes after post-processing mode  1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main
if size(need_split,1)>0
    clear post_segments  post_radii post_radiidx nodes;
    temp_count = 0;
    for n = 1:size(need_split,1)    
        Mask = L(:,:,:)== need_split(n);
        expname = 'results-r/10_exp'+string(n);
        need_splitn = need_split(n);
        [final_Groups,postV] = mLCuts(Mask,MODE,sizeLimit,distLimit,r,sigmaD,sigmaT,expname,L,need_splitn,showI,saveI);%, GT,tG,roi);
        for v = 1:size(unique(postV(:)),1)-1
            temp_count = temp_count+1;
            postV1(postV==v)=temp_count;
        end
        %pause;
        close all;
    end
end

%% Plot current reconstructed volume
no_need_split = [] ;
labels = sort(unique(L(:)));
num_e = size(unique(L(:)),1)-1;
count1 = 0;
Vwhole1 = zeros(size(L)); % reconstructed whole volume after first round
for n = 1:num_e
    temp_showV = zeros(size(L));
    temp_showV(L==labels(n+1)) = n;
    cc = regionprops3((temp_showV==n),'ConvexVolume');
    s = cc.ConvexVolume/sum(temp_showV==n,'all');
    if size(s,1) >1
        s = 0;
    end
    if sum(temp_showV==n,'all')< VOLUME_UPPERBOUND && s<= CONCAVITY_RATIO
       count1 = count1+1;
       Vwhole1(L==labels(n+1)) = count1;
       figure(2);isosurface(temp_showV,n-1);colormap colorcube;axis equal;title('Reconstructed volume after step 1')
       no_need_split = [no_need_split;labels(n+1)];
    end
end
for n = 1:temp_count
    temp_showV = zeros(size(L));
    temp_showV(postV1==n) = n;
    if sum(temp_showV==n,'all')> 0
       count1 = count1+1;
       Vwhole1(postV1==n) = count1;
       figure(2);isosurface(temp_showV,n-1);colormap colorcube;axis equal;axis off;
    end
end

%% Sequential step 2: find over-segmented masks
% The current selection includes the clusters that is smaller than a regular
% bacterial size; in addition, we will find their neighboring segments to
% decide if a connection is needed.
need2 = []; % indices of over-segmented clusters
needV = zeros(size(Vwhole1));
VOLUME_LOWERBOUND = 4/3*pi*4*4*20;
temp_count = 0;
for n = 1:count1
    temp_showV = zeros(size(L));
    currentseg = Vwhole1==n;
    if sum(currentseg(:))< VOLUME_LOWERBOUND
       temp_count = temp_count+1;
       temp_showV(Vwhole1==n) = n;
       needV(Vwhole1==n) = temp_count;
       need2 = [need2;n];
    end
end

% relabel and dilate needV to find its neighboring segments
[L2,num]  = bwlabeln(needV>0);
se = strel('cube',2);
LL = zeros(size(L2));
for r = 1:num
    LL(imdilate(L2 == r,se)) = r;
end

% multiple with input volume (Vwhole1) to find regions of Masks M for local merging;
% connected regions are masking with the same ID in M;
% this can make step 2 faster.
M = zeros(size(L2));
need_merge = []; %indices of all ROI's for merging
for m = 1:num
    temp_V = Vwhole1.*(LL==m);
    ids = unique(temp_V(:));
    for n = 1: size(unique(temp_V(:)),1)-1
        currentseg = Vwhole1 == ids(n+1); 
        if sum(currentseg(:)) <= 4/3*pi*4*16+pi*16*25 % 
            need_merge = [need_merge;ids(n+1)];
            M(Vwhole1 == ids(n+1))=m;          
        end
    end
end

%% Then solve oversegmentation: run m-LCuts mode 2
%%% Parameters and initializations %%%
showI = 0;
saveI = 0;
distLimit =25;
sigmaD = 5;
sigmaT = 0.1;
r = [1,7];
sizeLimit = [4/3*pi*4*4*40,4/3*pi*4*4*8];
MODE = 2; % there's a function in matlab called mode, so just change this paramter name to avoid conflict
postV2 = zeros(size(L));
clear final_Groups postV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_count = 0;
for i = 1:num
    %close all;
    expname = 'results-r/10_2exp'+string(i);
    curMask = (M==i).*Vwhole1; % ROI's in the current masking region
    temp_label = 0;
    Mask = zeros(size(curMask));
    ids = unique(curMask(:));
    for n = 1: size(unique(curMask(:)),1)-1
        temp_label = temp_label+1;
        Mask(curMask == ids(n+1))=temp_label; % relabel the cluster IDs in M==i
    end
    need_splitn = i;
    
    [final_Groups,postV] = mLCuts(Mask,MODE,sizeLimit,distLimit,r,sigmaD,sigmaT,expname,Vwhole1,need_splitn,showI,saveI);    
    for v = 1:size(unique(postV(:)),1)-1
        temp_count = temp_count+1;
        postV2(postV==v)=temp_count; % num of cells after post-processing
    end
end
%% Reconstruct the final volume
Vend = zeros(size(Vwhole1));% final reconstructed volume
count2 = 0;

for n = 1:count1
    temp_showV = zeros(size(L));
    currentseg = Vwhole1==n;
    if sum(currentseg(:))> VOLUME_LOWERBOUND && ~ismember(n,need_merge)
       count2 = count2+1;
       temp_showV(Vwhole1==n) = n;
       Vend(Vwhole1==n) = count2;
       figure(3);isosurface(temp_showV,n-1);colormap colorcube;axis equal;title('Reconstructed volume after step 2')
    end
end
for n = 1:temp_count
    temp_showV = zeros(size(L));
    temp_showV(postV2==n) = n;
    % to get the final reconstruction, one can also select to filter out
    % the too small noisy clutter
    if sum(temp_showV==n,'all')> 4/3*pi*4*4*4 
       count2 = count2+1;
       Vend(postV2==n) = count2;
       figure(3);isosurface(temp_showV,n-1);colormap colorcube;axis equal;axis off;
    end
end