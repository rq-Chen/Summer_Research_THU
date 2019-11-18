%% saveLabelData - save the data

% Ruiqi Chen, Aug. 15, 2019

clear;clc;close all;

%% For Rest State

clear;clc;close all;

%% Constants

DFPATH = '..\new\crq\rest\restRef\';
EEGFILE = 'crqRestRRData.mat';
load([DFPATH EEGFILE], 'eegdata');
load([DFPATH EEGFILE], 'eegdata', 'chanlocs');
% The file should contain the numeric matrix (channels * timepoints) of the
% pre-processed EEG data, entitled eegdata; and the cell array (channels * 1)
% of the name of each channel, entitled channels.
SAMPLERATE = 1000;  % (Hz)
WINLENRANGE = [10 20 45 100 200 700 2000];  % (Microsecond)
NWINLEN = size(WINLENRANGE, 2);
% load([DFPATH 'ClusteringData.mat'], 'ClusterData');
ClusterData = cell(4,NWINLEN);  % phi_0 = 0, 45, 90, 135

anaSig = nan(size(eegdata, 1), size(eegdata, 2));
hilLength = fix(SAMPLERATE * 0.5);  % number of points in 500ms
hilNWin = fix(size(eegdata, 2) / hilLength);
for i = 1:hilNWin - 1
    hilStart = hilLength * (i - 1) + 1;
    hilEnd = hilLength * i;
    anaSig(:, hilStart:hilEnd) = hilbert(eegdata(:, hilStart:hilEnd)')';
end
anaSig(:, hilEnd + 1:end) = hilbert(eegdata(:, hilEnd + 1:end)')';

for jLoop = 1:1
    for iLoop = 2:NWINLEN % should be 1:NWINLEN if newly created
        WINLEN = WINLENRANGE(iLoop);
        currStruct.Lag = fix(WINLEN / size(ClusterData, 1) * (jLoop - 1));
        currStruct.Winlen = WINLEN;
        currStruct.ClustData = cell(1, 12);
        tic;
        for K_CLUSTER = 3:5  % number of clusters

            %% Calculate the vector for clustering

            ptEachWin = fix(SAMPLERATE * WINLEN / 1000);  % number of timepoints within each window
            ptStart = 1 + currStruct.Lag;  % discard the data before it
            clNum = fix((size(eegdata, 2) - ptStart + 1) ...
                / ptEachWin);  % number of points for clustering


            clData = zeros(clNum, nchoosek(size(eegdata,1), 2));

            tic;
            for i = 1:clNum
                currentData = anaSig(:, (ptStart + (i - 1) * ptEachWin) : ...
                    (ptStart + i * ptEachWin - 1));  % get the data within this window
                currPhi = atan(imag(currentData) ./ real(currentData));
                currPhi = currPhi';
                loopI = 1;
                for icol = 1:size(currPhi, 2) - 1
                    for jcol = icol + 1:size(currPhi, 2)
                        clData(i, loopI) = PLI(currPhi(:, icol), currPhi(:, jcol));
                        loopI = loopI + 1;
                    end
                end                
%                 currentData = eegdata(:, (ptStart + (i - 1) * ptEachWin) : ...
%                     (ptStart + i * ptEachWin - 1));  % get the data within this window
%                 clData(i, :) = corrcoef(currentData');
%                 clData(i,:) = squareform(tril(cofMat,-1));
%                 % squareform doesn't check for symmetry, but does check for zero diag
            end
            toc;
            

            %% Clustering & Visualization

            [currClData.Idx, currClData.Centers] = kmeans(clData, K_CLUSTER, 'Replicates', 5);
            currClData.kVal = K_CLUSTER;
            currClData.clNum = size(clData, 1);

            % Silhouette
            currClData.Pt_Selected = 1:currClData.clNum;
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

            %% MDS visualization

%             disTri = pdist(clData);
%             currClData.MDS_Selected = mdscale(disTri, 2);

            currStruct.ClustData{1, K_CLUSTER} = currClData;
        end
        ClusterData{jLoop,iLoop} = currStruct;
        toc
    end
end
save([DFPATH 'PLIClusteringData.mat'], 'ClusterData');
save([DFPATH EEGFILE], 'eegdata', 'anaSig', 'chanlocs');

function pliVal = PLI(A, B)
    pliVal = mean(sign(A - B));
end
