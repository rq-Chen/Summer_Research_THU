%% cmpGFPfc.m - compare the outcome of two different clustering method

clear;clc;close all;

DFPATH = '../extra/';
CLUSTFILE = 'ClusteringData.mat';
GFPFILE = 'msGFP.mat';
OUTPUTFILE = 'cmpClustData.mat';
WINLENS = 2:5;  % change it if some winlens are not available
K_CLUSTERS = 3:5;
K_MAX = 6;
SAMPLERATE = 128;

load([DFPATH CLUSTFILE], 'ClusterData');
load([DFPATH GFPFILE], 'chanlocs', 'microstate');
GFPLabels = microstate.fit.labels;

%% Compare by counting

cntRes = nan(size(ClusterData, 2), max(K_CLUSTERS), K_MAX, K_MAX);

for i = WINLENS
    for j = K_CLUSTERS
        currClust = ClusterData{1, i};
        currPtWin = fix(currClust.Winlen * SAMPLERATE / 1000);
        Idx = currClust.ClustData{1, j}.Idx;
        perMat = zeros(max(GFPLabels), j);
        for k = 1:currPtWin * size(Idx, 1)
            iIdx = GFPLabels(k);
            jIdx = Idx(ceil(k / currPtWin));
            perMat(iIdx, jIdx) = perMat(iIdx, jIdx) + 1;
        end
        rFreq = sum(perMat, 2);
        rFreq = rFreq / sum(rFreq);
        cFreq = sum(perMat);
        cFreq = cFreq / sum(cFreq);
        estMat = sum(perMat, 'all') * (rFreq * cFreq);
        freqRatio = perMat ./ estMat;
        cntRes(i, j, 1:j, 1:max(GFPLabels)) = freqRatio';
    end
end

save([DFPATH OUTPUTFILE], 'cntRes');

%% Compare by correlation

% corrRes = nan(size(ClusterData, 2), max(K_CLUSTERS), K_MAX, K_MAX);
% 
% for i = 1:size(corrRes, 1)
%     for j = K_CLUSTERS
%
%     end
% end