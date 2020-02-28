%% IlluStTopo.m - illustrate the topology of the state center
%
% Ruiqi Chen, 2019

%% Parameters

DFPATH = '..\new\crq\rest\restref\';
CLFILE = 'PLIClusteringData.mat';
CHFILE = '..\avrRef\msGFP.mat';
WINLEN = [10 20 45 100 200 700 2000];
WINLENINDEX = 2:7;
K_CLUSTERS = 3:5;

%% Printing

eeglab;
figure;
load([DFPATH CLFILE], 'ClusterData');
load([DFPATH CHFILE], 'chanlocs');

for iTime = WINLENINDEX
    for iK = K_CLUSTERS
        Center = ClusterData{1, iTime}.ClustData{1, iK}.Centers;
        for i = 1:size(Center, 1)
            MDSVec = mdscale(logsig(-Center(i,:)), 1, 'Criterion', 'sstress');    
            topoplot(MDSVec, chanlocs);
            title(sprintf("State %d of %d", i, size(Center, 1)));
            saveas(gca, [DFPATH ...
                sprintf('PLIillustrate/%d/%d/%dTP.jpg', WINLEN(iTime), iK, i)]);
            clf;
        end
    end
end

close all;