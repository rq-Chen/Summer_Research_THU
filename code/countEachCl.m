%% countEachCl.m - count the propotion of each cluster
%
% Ruiqi Chen, 2019.10.7

DFPATH = '..\new\crq\rest\restref\';
CLFILE = 'ClusteringData.mat';
WINLENINDEX = 6;  % 700ms
K_CLUSTERS = 4;

load([DFPATH CLFILE], 'ClusterData');
ClusterData1 = ClusterData{1, WINLENINDEX}.ClustData{1, K_CLUSTERS}.Idx;
cnt = zeros(1, 4);
for i = 1:size(ClusterData1, 1)
    cnt(ClusterData1(i)) = cnt(ClusterData1(i)) + 1;
end
cnt = cnt / size(ClusterData1, 1)