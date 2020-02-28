%% clustByGFP.m - cluster the data with the segmentation by GFP
%
% INPUT:
%
% The EEGFILE should contain the numeric matrix (channels * timepoints) of the
% pre-processed EEG data, entitled eegdata; and the cell array (channels * 1)
% of the name of each channel, entitled channels.
%
% The GFPFILE should contain the chanlocs and microstate structure in the
% EEG structure from the EEGLAB, after performing GFP analysis with the
% microstates analysis plugin.

clear;clc;close all;

%% Constants

DFPATH = '..\new\hc\rest\avrRef\';
EEGFILE = 'hcRestAvrData.mat';
GFPFILE = 'msGFP.mat';
OUTPUTFILE = 'clustByGFP.mat';
load([DFPATH EEGFILE], 'eegdata');
load([DFPATH GFPFILE], 'microstate', 'chanlocs');
SAMPLERATE = 1000;  % (Hz)

%% Clustering

ClustData = cell(1, 12);
GFPLabels = microstate.fit.labels;

tic;
for K_CLUSTER = 3:5  % number of clusters
    
    % Calculate the vector for clustering
    
    clNum = 0;
    iStart = 1;
    iEnd = 2;
    clData = nan(size(GFPLabels, 2), nchoosek(size(eegdata,1), 2));
    
    while 1
        while iEnd <= size(GFPLabels, 2) ...
                && GFPLabels(iStart) == GFPLabels(iEnd)
            iEnd = iEnd + 1;
        end
        iEnd = iEnd - 1;
        clNum = clNum + 1;
        
        currentData = eegdata(:, iStart : iEnd);  % get the data within this window
        cofMat = corrcoef(currentData');
        clData(clNum,:) = squareform(tril(cofMat,-1));
        % squareform doesn't check for symmetry, but does check for zero diag
        
        if iEnd == size(GFPLabels, 2)
            break;
        end
        
        iStart = iEnd + 1;
        iEnd = iStart + 1;
    end
    clData = clData(1:clNum, :);
    
    % Clustering & Visualization
    
    [currClData.Idx, currClData.Centers] = kmeans(clData, K_CLUSTER, 'Replicates', 5);
    currClData.kVal = K_CLUSTER;
    currClData.clNum = clNum;
    
    % Silhouette
    currClData.Pt_Selected = 1:clNum;
    if currClData.clNum > 1600
        randPer = randperm(currClData.clNum);
        randPer = randPer(1:1600);
        currClData.Pt_Selected = randPer;
        clData = clData(randPer, :);
        idx = currClData.Idx(randPer,:);
    else
        idx = currClData.Idx;
    end
    currClData.Sil_Selected = silhouette(clData, idx);
    
    % MDS visualization
    
    %             disTri = pdist(clData);
    %             currClData.MDS_Selected = mdscale(disTri, 2);
    
    ClustData{1, K_CLUSTER} = currClData;
end
toc

save([DFPATH OUTPUTFILE], 'ClustData');

%% Visualiztion

for j = 3:5
    Pts = ClustData{1, j}.Pt_Selected;
    Sil = ClustData{1, j}.Sil_Selected;
    Labels = ClustData{1, j}.Idx;
    Labels = Labels(Pts);
    for l = 1:j
        sTitle = sprintf('state %d of %d, Sil val = %f', ...
            l, j, mean(Sil(Labels == l)));
        fprintf('%s\n', sTitle);
        
        heatMat = ClustData{1, j}.Centers;
        heatMat = squareform(heatMat(l,:));
        for x = 1 : size(heatMat, 1)
            heatMat(x,x) = 1;
        end
        figure;
        heatmap(heatMat, 'Colormap', parula);
        caxis([-1 1]);
        title(sTitle);
        colorbar off;
 
        folderName = [DFPATH 'illustrate/GFP/' num2str(j) '/'];
        mkdir(folderName);
        saveas(gca, [folderName sprintf('%d_%d.jpg', l, j)]);
        close all;
    end
end