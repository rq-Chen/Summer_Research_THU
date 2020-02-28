%% Illustrate - illustrate the distribution of microstates

% Ruiqi Chen, Aug. 14, 2019
% Illustrate the data by multidimensional scaling given the length of the
% interval.

%% Parameters
%
% eegdata: the numeric matrix (channels * timepoints) of the pre-processed 
% signal. If not specified, ClustSel reads it from the file entitled FILENAME.
%
% SAMPLERATE: in Hz, by default 200.
% 
% WINLEN: in microseconds, by default 200.
%

%% Function
function [] = Illustrate(eegdata, SAMPLERATE, WINLEN)

    tic;
    % Default value

    DFPATH = '..\\new\\ydy\\Language\\shantianfang\\';
    FILENAME = 'ydyStfData.mat';
    % The file contains the numeric matrix (channels * timepoints) of the
    % pre-processed EEG data, entitled eegdata.
    DFSAMRATE = 200;
    DFWINLEN = 200;
    
    % For phase, we only choose 0, 36, 72, 108, 144 degree (1/5 to 4/5)
    PTPHASE = [0 0.2 0.4 0.6 0.8];
    phaseN = size(PTPHASE, 2);
    
    if nargin < 1 || isempty(eegdata)
        load([DFPATH FILENAME], 'eegdata');
    end
    if nargin < 2 || isempty(SAMPLERATE)
        SAMPLERATE = DFSAMRATE;
    end
    if nargin < 3 || isempty(WINLEN)
        WINLEN = DFWINLEN;
    end
    ptEachWin = fix(SAMPLERATE * WINLEN / 1000);  % number of timepoints within each window

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
    
    % Transforming
    for i = 1:phaseN
        for j = 1:clNum(i)
            currentData = eegdata(:, (iLag(i) + 1 + (j - 1) * ptEachWin) : ...
                (iLag(i) + j * ptEachWin));  % get the data within this window
            cofMat = corrcoef(currentData');
            clData(i,j,:) = squareform(tril(cofMat, -1));
        end
    end

    % MDS visualization
    tempP = gcp('nocreate');
    if ~isempty(tempP)
        delete(tempP);
    end
    parpool(min([phaseN 4]));
    parfor i = 1:phaseN
        currData = squeeze(clData(i,:,:));
        disTri = pdist(currData(1:clNum(i),:));
        MDSMat = mdscale(disTri, 2);

        % create subfolder first
        mkdir([DFPATH sprintf('%d', WINLEN)]);
        
        figure;
        scatter(MDSMat(:,1), MDSMat(:,2));

        title(sprintf("winlen = %d, lag = %d", WINLEN, iLag(i)));
        savefig([DFPATH sprintf('%d\\lag%d.fig', WINLEN, iLag(i))]);
    end
    
    toc
end