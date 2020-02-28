%% parClustSel - the paralleled ClustSel function

% Ruiqi Chen, Aug. 12, 2019
% Input the EEG data and a specific window length (e.g. 200ms), calculate
% the best k (number of clusters) and lag (or the beginning phase), then
% illustrate the best result (optional).
%
% Caution: this program runs very slow, because although the calculation of
% correlation coefficient matrix and k-means clustering are both fast, it
% takes nearly 50s to compute the silhouette value for a specific solution.
% Currently, we only parallelize the outer loop (for phase) and use the
% originally evalclusters() function. Since we only have 4 cores but 5
% phases, it is unclear whether the process can be accelerated
% significantly. The next step is to break out the evalclusters() function
% and parallelize at the level of silhouette value computation.

%% Parameters
%
% eegdata: the numeric matrix (channels * timepoints) of the pre-processed 
% signal. If not specified, ClustSel reads it from the file entitled FILENAME.
%
% SAMPLERATE: in Hz, by default 200.
% 
% WINLEN: in microseconds, by default 200.
%
% DEPICT: by default 0 (no visualization). If DEPICT == 1, Silhouette plot;
% if DEPICT == 2, MDS visualization; if DEPICT == 3, both.

%% Return Values
%
% optK: the optimal number of clusters
%
% optLag: the optimal lag (counted by timepoints), namely the number of
% timepoints discarded at the front of the data
%
% optIdx: the correspondent clustering solution

%% Function
function [optK, optLag, optIdx, optData] = parClustSel(eegdata, SAMPLERATE, WINLEN, DEPICT)

    % Default value

    DFPATH = '..\\new\\ydy\\Language\\shantianfang\\';
    FILENAME = 'ydyStfData.mat';
    % The file contains the numeric matrix (channels * timepoints) of the
    % pre-processed EEG data, entitled eegdata.
    DFSAMRATE = 200;
    DFWINLEN = 200;
    DFDEPICT = 0;
    
    % For phase, we only choose 0, 36, 72, 108, 144 degree (1/5 to 4/5)
    PTPHASE = [0 0.2 0.4 0.6 0.8];
    phaseN = size(PTPHASE, 2);
    KLIST = 3:8;  % list of number-of-clusters in interest
    
    if nargin < 1 || isempty(eegdata)
        load([DFPATH FILENAME], 'eegdata');
    end
    if nargin < 2 || isempty(SAMPLERATE)
        SAMPLERATE = DFSAMRATE;
    end
    if nargin < 3 || isempty(WINLEN)
        WINLEN = DFWINLEN;
    end
    if nargin < 4 || isempty(DEPICT)
        DEPICT = DFDEPICT;
    end
    ptEachWin = fix(SAMPLERATE * WINLEN / 1000);  % number of timepoints within each window

    
    % Calculate the best phase and k
    
    % definition of variables
    clNum = nan(phaseN, 1);
    % number of points for clustering (for a specific phase, same below)
    iLag = clNum;  % lag (counted in timepoints)
    for i = 1:phaseN
        iLag(i) = fix(PTPHASE(i) * ptEachWin);
        clNum(i) = fix((size(eegdata, 2) - iLag(i)) / ptEachWin);
    end
    maxClNum = max(clNum);
    clData = nan(phaseN, maxClNum,...
        nchoosek(size(eegdata,1), 2));  % Data
    crtVal = nan(phaseN, 1);  % the opt silhouette value
    crtK = crtVal;  % the opt number of clusters
    crtIdx = nan(phaseN, maxClNum);  % the opt solution
    
    
    % Transforming
    for i = 1:phaseN
        for j = 1:clNum(i)
            currentData = eegdata(:, (iLag(i) + 1 + (j - 1) * ptEachWin) : ...
                (iLag(i) + j * ptEachWin));  % get the data within this window
            cofMat = corrcoef(currentData');
            clData(i,j,:) = squareform(tril(cofMat, -1));
        end
    end
    
    % Clustering & Selection
    tic
    tempP = gcp('nocreate');
    if ~isempty(tempP)
        delete(tempP);
    end
    parpool(min([phaseN 4]));
    parfor i = 1:phaseN
        currData = squeeze(clData(i,:,:));
        eva = evalclusters(currData(1:clNum(i),:), 'kmeans', ...
            'silhouette', 'klist', KLIST);
        crtK(i) = eva.OptimalK;
        tmp = nan(maxClNum, 1);
        tmp(1:clNum(i)) = eva.OptimalY;
        crtIdx(i, :) = tmp;
        crtVal(i) = eva.CriterionValues(find(KLIST == eva.OptimalK, 1));
    end
    ParTime = toc
    [optVal, optI] = max(crtVal);
    optK = crtK(optI);
    optLag = iLag(optI);
    optIdx = squeeze(crtIdx(optI, 1:clNum(optI)));
    optData = squeeze(clData(optI,1:clNum(optI), :));

    % Plot the outcome
    tic
    if DEPICT > 0

        % create subfolder first
        mkdir([DFPATH sprintf('%d', WINLEN)]);
        
        figure;
        silhouette(optData, optIdx);

        hdl = gca;
        hdl.Children.EdgeColor = [.8 .8 1];
        xlabel 'Silhouette Value'
        ylabel 'Cluster'
        title(sprintf("k = %d, lag = %d, silhouette val = %f",...
            optK, optLag, optVal));
        savefig([DFPATH sprintf('%d\\optSil.fig', WINLEN)]);
    end

    % MDS visualization
    if DEPICT > 1
        
        disTri = pdist(optData);
        MDSMat = mdscale(disTri, 2);

        figure;
        gscatter(MDSMat(:,1), MDSMat(:,2), optIdx);
        
        title(sprintf("k = %d, lag = %d, silhouette val = %f",...
            optK, optLag, optVal));
        savefig([DFPATH sprintf('%d\\optMds.fig', WINLEN)]);
    end
    
    % heatmap of the center
    if DEPICT > 2
        for i = 1:optK
            currState = mean(optData(optIdx == i, :));
            currState = squareform(currState);
            for j = 1:size(currState,1)
                currState(j,j) = 1;
            end
            heatmap(currState, 'Colormap', parula);
            caxis([-1 1]);
            savefig([DFPATH sprintf('%d\\State%d.fig',...
                WINLEN, i)]);
            close all;
        end
    end
    toc
end