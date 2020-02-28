%% IlluSingle.m - show a single connectivity matrix

% Ruiqi Chen, Aug 19, 2019

clear; close;

%% Constants

DFPATH = '..\new\crq\Language\shantianfang\avrRef\';
EEGFILE = 'crqAvrData.mat';
% The file should contain the numeric matrix (channels * timepoints) of the
% pre-processed EEG data, entitled eegdata; and the cell array (channels * 1)
% of the name of each channel, entitled channels. (chanlocs is also OK)
CLFILE = 'ClusteringData.mat';
% The file should contain the struct cell ClusteringData

SAMPLERATE = 1000;  % (Hz)
WINLENINDEX = 2:5;  % index
K_CLUSTERS = 3:5;

%% Calculate the vector for clustering

load([DFPATH EEGFILE], 'eegdata');
load([DFPATH CLFILE], 'ClusterData');
figure;

for i = WINLENINDEX  % 10 - 200ms
    WINLEN = ClusterData{1, i}.Winlen;
    LAG = ClusterData{1, i}.Lag;        
    for j = K_CLUSTERS
        currClust = ClusterData{1, i}.ClustData{1, j};

        ptEachWin = fix(SAMPLERATE * WINLEN / 1000);  % number of timepoints within each window
        clNum = fix((size(eegdata, 2) - LAG) ...
            / ptEachWin);  % number of points for clustering

        clData = zeros(clNum, nchoosek(size(eegdata, 1), 2));

        %% Calculate the vector for clustering        
        tic
        for k = 1:clNum
            currentData = eegdata(:, (LAG + 1 + (k - 1) * ptEachWin) : ...
                (LAG + k * ptEachWin));  % get the data within this window
            cofMat = corrcoef(currentData');
            clData(k,:) = squareform(tril(cofMat,-1));
        end
        
        SampleData = zeros(j, 10, size(clData, 2));
        Idx = currClust.Idx;
        for k = 1:j
            tmp = find(Idx == k);
            randp = randperm(size(tmp, 1));

            folderName = [DFPATH 'illustrate/' num2str(WINLEN)...
                '/' num2str(j) '/' num2str(k) '/'];
            mkdir(folderName);
            
            histogram(currClust.corrData(k, tmp));
            sTitle = sprintf('winlen = %d, state %d of %d, mean = %f', ...
                WINLEN, k, j, mean(currClust.corrData(k, tmp)));
            title(sTitle);
            saveas(gca, [folderName sprintf('%dCorr%d_%d.jpg',...
                WINLEN, k, j)]);
            
            SampleData(k, :, :) = clData(tmp(randp(1:10)),:);
            for l = 1:10
                tmp = squeeze(SampleData(k, l, :));
                tmp = squareform(tmp);
                for h = 1:size(tmp, 1)
                    tmp(h,h) = 1;
                end
                heatmap(tmp, 'Colormap', parula);
                caxis([-1 1]); colorbar off;
                saveas(gca, [folderName sprintf('%dCorr%d_%d_%d.jpg',...
                    WINLEN, k, j, l)]);
            end
        end
        
        ClusterData{1,i}.ClustData{1,j}.Sample = SampleData;
    end
end

save([DFPATH CLFILE], 'ClusterData');
close all;