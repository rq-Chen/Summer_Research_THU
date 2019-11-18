%% corrLanClust.m - allign the clustering and linguistic label

% Ruiqi Chen, Aug. 21, 2019
% Calculate the frequency of timepoints within each word interval being clustered
% to different groups. Further processing may include normalizing the
% frequency by the independence-assumed frequency.

clear;close all;clc;

%% Parameters

DFPATH = '../new/ydy/Language/6_13/avrRef/';
CLFILE = 'remove/reClusteringData.mat';
LAGINDEX = 1;  % first row, lag = 0
WINLENINDEX = 2;  % which column
CL_OI = 1;  % the cluster of interest
K_CLUSTERS = 4;
LNFILE = 'ydy13Labels.mat';
EEGFILE = 'ydy13AvrData.mat';
OUTPUTFILE = 'remove/ydy13LanLabels.mat';
SAMPLERATE = 1000;
MODE = 1;

% CLFILE contains the variable ClusterData (4 * 7 struct cell, but only the
% first row (where lag = 0) will be used). LNFILE contains the variable
% MusicLabel (N * 4 string cell, onset offset category text, and will be
% extended by K_CLUSTERS more columns containing the counting number we
% want). EEGFILE contains the eegdata (n_channels * n_timepoints double) and
% channels (n_channels * 1 cell). All the data must be alligned beforehand.

% MODE = 1: select by relative frequency; MODE = 2: select by absolute
% frequency; MODE = 3: select by correlation.

%% Count the frequency

load([DFPATH CLFILE]);
load([DFPATH LNFILE]);
load([DFPATH EEGFILE]);

LanClustLabel = cell(size(MusicLabel, 1), size(MusicLabel, 2) + K_CLUSTERS);
LanClustLabel(:, 1:size(MusicLabel, 2)) = MusicLabel;
LanClustLabel(:, size(MusicLabel, 2) + 1 : end) = ...
    num2cell(zeros(size(MusicLabel, 1), K_CLUSTERS));

Clust = ClusterData{LAGINDEX, WINLENINDEX}.ClustData{1, K_CLUSTERS};
Centers = Clust.Centers;
Idx = Clust.Idx;
Corr = Clust.corrData;
CorrSort = sort(Corr, 2);
CorrThres = CorrSort(:, fix(size(Corr, 2) * 0.99));
TpIdx = nan(1, size(eegdata, 2));
Lag = ClusterData{LAGINDEX, WINLENINDEX}.Lag;
Winlen = ClusterData{LAGINDEX, WINLENINDEX}.Winlen;

for i = 1:size(Idx, 1)
    if MODE ~= 3
        tmp = Idx(i);
    else
        tmp = find((Corr(:, i) >= CorrThres), 1);
    end
    
    if tmp
        TpIdx(Lag + 1 + (i - 1) * Winlen : Lag + i * Winlen) = ...
            ones(1, Winlen) * tmp;
    end
end
TpIdx = categorical(TpIdx);

for i = 1:size(LanClustLabel, 1)
    startIdx = ceil(LanClustLabel{i, 1} * SAMPLERATE);
    endIdx = ceil(LanClustLabel{i, 2} * SAMPLERATE);
    if endIdx > size(eegdata, 2)
        LanClustLabel = LanClustLabel(1:i - 1, :);
        break;
    end
    LanClustLabel(i, size(MusicLabel, 2) + 1 : end) = ...
        num2cell(countcats(TpIdx(startIdx:endIdx)));
end

nCount = cell2mat(LanClustLabel(:, size(MusicLabel, 2) + 1 : end));
if MODE == 1 || MODE == 3
    rFreq = sum(nCount, 2);
    rFreq = rFreq / sum(rFreq);
    cFreq = countcats(TpIdx);
    cFreq = cFreq / sum(cFreq, 2);
    FreqEst = rFreq * cFreq;
    FreqEst = FreqEst * sum(nCount, 'all');
    FrqRatio = nCount ./ FreqEst;
    [~, MaxIdx] = max(FrqRatio, [], 2);
    MaxIdx(isnan(FrqRatio(:,1))) = nan(sum(isnan(FrqRatio(:, 1))), 1);
end

if MODE == 2
    [~, MaxIdx] = max(nCount, [], 2);
end

% if MODE == 3
%     MaxIdx = cell2mat(LanClustLabel(:, size(MusicLabel, 2) + 1)) > 0;
% end


wordOI = LanClustLabel(MaxIdx == CL_OI, :);

allWords = LanClustLabel(:, 3);
for i = 1:size(allWords, 1)
    if MaxIdx(i) == CL_OI
        allWords{i, 1} = join(['**' allWords{i, 1} '**'], '');
    end
end
allWords = join([allWords{:}]);

save([DFPATH OUTPUTFILE], 'LanClustLabel', 'Centers', 'Lag', 'Winlen',...
    'wordOI', 'allWords');



