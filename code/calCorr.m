%% calCorr.m - calculate the correlation curve

% Ruiqi Chen, Aug 19, 2019
% Calculate the correlation between the activity within each time window
% and the center of the cluster of interest.

%% Constants

DFPATH = '..\new\wrx\shantianfang\Cz\';
EEGFILE = 'wrxStfCzData.mat';
% The file should contain the numeric matrix (channels * timepoints) of the
% pre-processed EEG data, entitled eegdata; and the cell array (channels * 1)
% of the name of each channel, entitled channels.
CLFILE = 'ClusteringData.mat';
% The file should contain the struct cell ClusteringData

SAMPLERATE = 1000;  % (Hz)
WINLENINDEX = 2:5;  % index
K_CLUSTERS = 3:5;

%% Calculate the vector for clustering

load([DFPATH EEGFILE], 'eegdata');
load([DFPATH CLFILE], 'ClusterData');

for i = WINLENINDEX  % 10 - 200ms
    WINLEN = ClusterData{1, i}.Winlen;
    LAG = ClusterData{1, i}.Lag;        
    for j = K_CLUSTERS
        tmp = ClusterData{1, i}.ClustData{1, j}.Centers;
        Centers = nan(size(tmp, 1) + 1, size(tmp, 2));
        Centers(1:end - 1, :) = tmp;

        ptEachWin = fix(SAMPLERATE * WINLEN / 1000);  % number of timepoints within each window
        clNum = fix((size(eegdata, 2) - LAG) ...
            / ptEachWin);  % number of points for clustering

        corrData = zeros(j, clNum);

        %% Calculate the vector for clustering        
        tic
        for k = 1:clNum
            currentData = eegdata(:, (LAG + 1 + (k - 1) * ptEachWin) : ...
                (LAG + k * ptEachWin));  % get the data within this window
            cofMat = corrcoef(currentData');
            tmp = squareform(tril(cofMat,-1));
            Centers(end, :) = tmp;
            cofMat = corrcoef(Centers');
            corrData(:, k) = cofMat(end, 1:end - 1);
        end
        ClusterData{1,i}.ClustData{1,j}.corrData = corrData;
    end
end

save([DFPATH CLFILE], 'ClusterData');