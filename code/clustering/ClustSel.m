%% ClustSel - cluster and find the best k given winlen and samprate

% Ruiqi Chen, Aug. 12, 2019
% Input the EEG data and a specific window length (e.g. 200ms), calculate
% the best k (number of clusters) and lag (or the beginning phase), then
% illustrate the best result.

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
%
% optData: the data computed under optK and optLag

%% Function
function [optK, optLag, optIdx, optData] = ClustSel(eegdata, SAMPLERATE, WINLEN, DEPICT)

    tic
    % Default value

    FILENAME = '..\new\ydy\Language\shantianfang\ydyStfData.mat';
    % The file contains the numeric matrix (channels * timepoints) of the
    % pre-processed EEG data, entitled eegdata.
    DFSAMRATE = 200;
    DFWINLEN = 200;
    DFDEPICT = 0;
    KILST = 3:12;
    crtVal = NaN;
    
    if nargin < 1 || size(eegdata, 1) == 0
        load(FILENAME, 'eegdata');
    end
    if nargin < 2
        SAMPLERATE = DFSAMRATE;
    end
    if nargin < 3
        WINLEN = DFWINLEN;
    end
    if nargin < 4
        DEPICT = DFDEPICT;
    end
    ptEachWin = fix(SAMPLERATE * WINLEN / 1000);  % number of timepoints within each window

    
    % Calculate the best phase and k
    
    % For phase, we only choose 0, 36, 72, 108, 144 degree (1/5 to 4/5)
    PTPHASE = [0 0.2 0.4 0.6 0.8];

    for iLag = fix(PTPHASE * ptEachWin)
        
        clNum = fix((size(eegdata, 2) - iLag) / ptEachWin);  % number of points for clustering
        clData = zeros(clNum, nchoosek(size(eegdata,1), 2));  % Data

        % Transforming
        for i = 1:clNum
            currentData = eegdata(:, (iLag + 1 + (i - 1) * ptEachWin) : ...
                (iLag + i * ptEachWin));  % get the data within this window
            cofMat = corrcoef(currentData');
            clData(i,:) = squareform(tril(cofMat, -1));
        end


        % Clustering & Selection
        eva = evalclusters(clData, 'kmeans', 'silhouette', 'klist', KILST);
        if isnan(crtVal) || eva.CriterionValues(eva.OptimalK) > crtVal
            optData = clData;
            optK = eva.OptimalK;
            optIdx = eva.OptimalY;
            optLag = iLag;
            crtVal = eva.CriterionValues(find(KILST == optK, 1));
        end

    end

    
    % Plot the outcome
    if DEPICT > 0

        % create subfolder first
        mkdir(sprintf('..\\new\\ydy\\Language\\shantianfang\\%d', ...
            WINLEN));
        
        figure;
        silhouette(optData, optIdx);

        hdl = gca;
        hdl.Children.EdgeColor = [.8 .8 1];
        xlabel 'Silhouette Value'
        ylabel 'Cluster'
        title(sprintf("k = %d, lag = %d, silhouette val = %f",...
            optK, optLag, crtVal));
        savefig(['..\\new\\ydy\\Language\\shantianfang\\' ...
            sprintf('%d\\optSil.fig', WINLEN)]);
    end

    % MDS visualization
    if DEPICT == 2
        
        disTri = pdist(optData);
        MDSMat = mdscale(disTri, 2);

        figure;
        gscatter(MDSMat(:,1),MDSMat(:,2), optIdx);
        
        title(sprintf("k = %d, lag = %d, silhouette val = %f",...
            optK, optLag, crtVal));
        savefig(['..\\new\\ydy\\Language\\shantianfang\\' ...
            sprintf('%d\\optMds.fig', WINLEN)]);
    end
    toc
end