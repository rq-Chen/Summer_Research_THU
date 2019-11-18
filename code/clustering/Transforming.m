%% Transforming - transform the EEG signal for clustering

% Ruiqi Chen, Aug. 12, 2019

clear;clc;close all;

%% Constants

EEGFILE = '..\new\ydy\Language\shantianfang\ydyStfData.mat';
% The file should contain the numeric matrix (channels * timepoints) of the
% pre-processed EEG data, entitled eegdata; and the cell array (channels * 1)
% of the name of each channel, entitled channels.
SAMPLERATE = 1000;  % (Hz)
WINLEN = 0.1;  % (second)
K_CLUSTER = 3;  % number of clusters

%% Calculate the vector for clustering

% load(EEGFILE);
ptEachWin = fix(SAMPLERATE * WINLEN);  % number of timepoints within each window
ptStart = 1;  % discard the data before it
clNum = fix((size(eegdata, 2) - ptStart + 1) ...
    / ptEachWin);  % number of points for clustering


clData = zeros(clNum, nchoosek(size(eegdata,1), 2));

tic
for i = 1:clNum
    currentData = eegdata(:, (ptStart + (i - 1) * ptEachWin) : ...
        (ptStart + i * ptEachWin - 1));  % get the data within this window
    cofMat = corrcoef(currentData');
    clData(i,:) = squareform(tril(cofMat,-1));
    % squareform doesn't check for symmetry, but does check for zero diag
    
%     if mod(i,200) == 0  % illustrate a sample result
%         heatmap(cofMat, 'Colormap', parula);
%         caxis([-1 1]); colorbar off;
%         DFPATH = '../new/hc/rest/avrRef/illustrate/GFPConcate/MS2MS5/';
%         title(sprintf('HC connectivity of GFP-based concatenated data, No. %d', ...
%             i / 100));
%         saveas(gcf, [DFPATH sprintf('%d.jpg', i / 100)]);
%     end
end
toc


%% Clustering & Visualization

[idx, Centers] = kmeans(clData, K_CLUSTER, 'Replicates', 5);

% Plot the outcome
tic;
close all;figure;
[sil, hdl] = silhouette(clData, idx);
ElapseTime = toc
sil = mean(sil);

hdl = gca;
hdl.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value'
ylabel 'Cluster'
title(sprintf("k = %d, average sil val = %f", K_CLUSTER, sil));


%% MDS visualization

tic
disTri = pdist(clData);
opts = statset('Display', 'final');
MDSMat = mdscale(disTri, 2, 'Options', opts);
ElapseTime = toc;

figure;
gscatter(MDSMat(:,1),MDSMat(:,2),idx);
title(sprintf("k = %d, average sil val = %f", K_CLUSTER, sil));
