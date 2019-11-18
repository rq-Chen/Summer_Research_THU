%% reCluster - remove a component from the data and re-cluster

% Ruiqi Chen, Aug. 26, 2019

%% Constants

DFPATH = '..\new\ydy\Language\6_13\avrRef\';
EEGFILE = 'ydy13AvrData.mat';
% The file should contain the numeric matrix (channels * timepoints) of the
% pre-processed EEG data, entitled eegdata; and the cell array (channels * 1)
% of the name of each channel, entitled channels.
CLFILE = 'ClusteringData.mat';
mkdir([DFPATH 'remove\']);
OUTPUTFILE = 'remove\reClusteringData.mat';
% The file should contain the struct cell ClusteringData

SAMPLERATE = 1000;  % (Hz)
WINLENINDEX = 1:5;  % index of the winlen
K_CLUSTERS = 4;
CLREMOVED = [4 2 1 2 3];  % No. of state to remove from data

%% Calculate the vector for clustering

load([DFPATH EEGFILE], 'eegdata');
load([DFPATH CLFILE], 'ClusterData');
reClusterData = cell(size(ClusterData, 1), size(ClusterData, 2));

for i = WINLENINDEX

    WINLEN = ClusterData{1, i}.Winlen;
    LAG = ClusterData{1, i}.Lag;     
    Centers = ClusterData{1, i}.ClustData{1, K_CLUSTERS}.Centers;
    ClRmv = Centers(CLREMOVED(i), :);

    currStruct.Lag = LAG;
    currStruct.Winlen = WINLEN;
    currStruct.ClustData = cell(1, 12);

    ptEachWin = fix(SAMPLERATE * WINLEN / 1000);  % number of timepoints within each window
    clNum = fix((size(eegdata, 2) - LAG) ...
        / ptEachWin);  % number of points for clustering
    mfClData = zeros(clNum, nchoosek(size(eegdata, 1), 2));

    %% Calculate the vector for clustering        

    for k = 1:clNum
        currentData = eegdata(:, (LAG + 1 + (k - 1) * ptEachWin) : ...
            (LAG + k * ptEachWin));  % get the data within this window
        cofMat = corrcoef(currentData');
        tmp = squareform(tril(cofMat,-1));
        mfClData(k, :) = tmp - ClRmv * (ClRmv * tmp') / (ClRmv * ClRmv');
    end

    %% Clustering & Visualization

    for j = 3:5

        [currClData.Idx, currClData.Centers] = kmeans(mfClData, j, 'Replicates', 5);
        currClData.kVal = j;
        currClData.clNum = size(mfClData, 1);

        % Silhouette
        currClData.Pt_Selected = 1:currClData.clNum;
        if currClData.clNum > 1600
            randPer = randperm(currClData.clNum);
            randPer = randPer(1:1600);
            currClData.Pt_Selected = randPer;
            ClData = mfClData(randPer, :);
            idx = currClData.Idx(randPer,:);
        else
            idx = currClData.Idx;
            ClData = mfClData;
        end
        currClData.Sil_Selected = silhouette(ClData, idx);

        %% MDS visualization

        disTri = pdist(ClData);
        currClData.MDS_Selected = mdscale(disTri, 2);

        currStruct.ClustData{1, j} = currClData;
    end
    reClusterData{1, i} = currStruct;
end

ClusterData = reClusterData;
save([DFPATH OUTPUTFILE], 'ClusterData');