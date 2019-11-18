%% IlluTopoChange.m - show a movie of the topology
%
% Ruiqi Chen, 2019.9.22
%
% For mode 1, please open eeglab before running this script.
% For mode 3, you need to change the printing function at the last of the
% script.

close all;clc;

%% Arguments

MODE = 3;
% select 1 for voltage value, 2 for functional connectivity patern,
% 3 for functional connectivity correlation

DFPATH = '../new/crq/music/';
EEGFILE = 'crqMusicOData.mat';
CLFILE = 'ClusteringData.mat';
CHFILE = '../rest/avrRef/msGFP.mat';  % file containing the chanlocs and microstate struct
OUTPUTFILE = 'illustrate/fcChange_Norm_200_in_2000.gif';

INTERVAL = 200;  % in microseconds
SAMPLERATE = 1000;  % in Hz
FPS = 3;  % frames per second
NTIME = 30;  % length of the output gif, in seconds

% for MODE == 1
SCALE = [-20 20];  % in uV
% for MODE == 2
RADIUS = 0.5150;
PROTHRES = 0.05;  % illustrate strongest connection
% for MODE == 3
WINLEN = [10 20 45 100 200 700 2000];
WINLENINDEX = 7;  % compare with the center of this winlen
K_CLUSTERS = 4;

load([DFPATH CHFILE], 'chanlocs', 'microstate');
allTheta = cell2mat({chanlocs(:).theta}) + 90;
allR = cell2mat({chanlocs(:).radius});
allX = allR .* cosd(allTheta);
allY = allR .* sind(allTheta);
circleX = RADIUS .* cosd(0:5:360);
circleY = RADIUS .* sind(0:5:360);

%% Generate pictures

GFPPeaks = sort(microstate.GFPpeakidx);

if MODE == 1
    
    load([DFPATH EEGFILE], 'eegdata');
    nFrames = FPS * NTIME;
    ptPerFrame = fix(SAMPLERATE * INTERVAL / 1000);
    selRange = nFrames * ptPerFrame + 1 : ...
        size(eegdata, 2) - nFrames * ptPerFrame - 1;
    ptStart = selRange(ceil(rand(1) * size(selRange, 2)));
    
    img = cell(1, nFrames);
    fig = figure;
    for i = 1:nFrames
        currPt = ptStart + (i - 1) * ptPerFrame;
        currVol = eegdata(:, currPt);
        
        isPeak = 0;
        for j = 1:size(GFPPeaks, 2)
            if GFPPeaks(j) < currPt - ptPerFrame / 2
                continue;
            end
            if GFPPeaks(j) < currPt + ptPerFrame / 2
                isPeak = 1;
            end
            break;
        end
        
        fig = clf;
        topoplot(currVol, chanlocs, 'maplimits', SCALE);
        if isPeak
            timeInt = 1;
            title('GFP peak');
        else
            timeInt = 1 / FPS;
        end
        img{1, i} = frame2im(getframe(fig));
        [A,map] = rgb2ind(img{1, i},256);
        
        if i == 1
            imwrite(A,map,[DFPATH OUTPUTFILE],'gif','LoopCount',Inf,...
                'DelayTime', timeInt);
        else
            imwrite(A,map,[DFPATH OUTPUTFILE],'gif','WriteMode','append',...
                'DelayTime', timeInt);
        end
    end
    close;
    
elseif MODE == 2
    
    load([DFPATH EEGFILE], 'eegdata');
    nFrames = FPS * NTIME;
    ptPerFrame = fix(SAMPLERATE * INTERVAL / 1000);
    selRange = nFrames * ptPerFrame + 1 : ...
        size(eegdata, 2) - nFrames * ptPerFrame - 1;
    ptStart = selRange(ceil(rand(1) * size(selRange, 2)));
    
    img = cell(1, nFrames);
    fig = figure;
    for i = 1:nFrames
        currStart = ptStart + (i - 1) * ptPerFrame;
        currEnd = currStart + ptPerFrame - 1;
        currData = eegdata(:, currStart:currEnd);
        cofMat = corrcoef(currData');
        
        isPeak = 0;
        for j = 1:size(GFPPeaks, 2)
            if GFPPeaks(j) < currStart
                continue;
            end
            if GFPPeaks(j) <= currEnd
                isPeak = 1;
            end
            break;
        end
        
        fig = clf;
        % %         imagesc(cofMat,[-1 1]);
        % %         caxis([-1 1]);
        %         topoplot(mdsVec, chanlocs);
        
        tmp = sort(squareform(tril(cofMat, -1)), 'descend');
        THRES = tmp(ceil(PROTHRES * size(tmp, 2)));
        
        hold on;
        for ll = 1:size(cofMat, 1)
            plot(circleX, circleY, 'k');
            plot(allX, allY, '.k', 'MarkerSize', 12);
        end
        for ll = 1:size(cofMat, 1)
            for lll = ll + 1 : size(cofMat, 1)
                if cofMat(ll, lll) > THRES
                    tmp = abs(cofMat(ll, lll));
                    if tmp < 0.85
                        colorTmp = [1 1 0];  % yellow
                    elseif tmp < 0.95
                        colorTmp = [0.9290 0.6940 0.1250];  % matlab yellow
                    else
                        colorTmp = [1 0 0];  % red
                    end
                    plot(allX([ll, lll]), allY([ll, lll]), ...
                        'Color', colorTmp, 'LineWidth', ...
                        2 * tmp);
                end
            end
        end
        axis equal; axis off;
        
        if isPeak
            timeInt = 1;
            title('GFP peak');
        else
            timeInt = 1 / FPS;
        end
        img{1, i} = frame2im(getframe(fig));
        [A,map] = rgb2ind(img{1, i},256);
        
        if i == 1
            imwrite(A,map,[DFPATH OUTPUTFILE],'gif','LoopCount',Inf,...
                'DelayTime', timeInt);
        else
            imwrite(A,map,[DFPATH OUTPUTFILE],'gif','WriteMode','append',...
                'DelayTime', timeInt);
        end
    end
    close;
    
else
    
    load([DFPATH CLFILE], 'ClusterData');
    load([DFPATH EEGFILE], 'eegdata');
    Labels = ClusterData{1, WINLENINDEX}.ClustData{1, K_CLUSTERS}.Idx;
    tmp = ClusterData{1, WINLENINDEX}.ClustData{1, K_CLUSTERS}.Centers;
    [CentersCor, csIndex] = sort(tmp, 2, 'descend');
    csIndex = csIndex(:, 1:ceil(size(tmp, 2) * PROTHRES));    
    CentersCor = CentersCor(:, 1:ceil(size(tmp, 2) * PROTHRES));
%     CentersCor = tmp(:, csIndex(1,:)); 

    
    nFrames = FPS * NTIME;
    ptPerFrame = fix(SAMPLERATE * INTERVAL / 1000);
    ptStart = 1;
    
    allCor = zeros(K_CLUSTERS, nFrames);
    for i = 1:nFrames
        currStart = ptStart + (i - 1) * ptPerFrame;
        currEnd = currStart + ptPerFrame - 1;
        currData = eegdata(:, currStart:currEnd);
        for j = 1:K_CLUSTERS
            tmp = squareform(tril(corrcoef(currData'),-1));
            tmp = corrcoef(CentersCor(j,:), tmp(csIndex(j,:)));
%             tmp = corrcoef(CentersCor(j,:), tmp(csIndex(1,:)));
            allCor(j, i) = tmp(1, 2);
        end
%         
%         allCor(:, i) = allCor(:, i) - mean(allCor(:, i));
%         
        
    end
    
    img = cell(1, nFrames);
    fig = figure;
    for i = 1:nFrames
        isPeak = 0;
        for j = 1:size(GFPPeaks, 2)
            currStart = ptStart + (i - 1) * ptPerFrame;
            currEnd = currStart + ptPerFrame - 1;
            if GFPPeaks(j) < currStart
                continue;
            end
            if GFPPeaks(j) <= currEnd
                isPeak = 1;
            end
            break;
        end
        
        fig = clf;
        xBegin = max(1, i - 20);
        xEnd = max(20, i);
        plot(xBegin:i, allCor(:, xBegin:i));
        axis([xBegin xEnd -0.5 0.5]);
        legend;
        
        MSNum = ceil(i * INTERVAL / WINLEN(WINLENINDEX));
        title(sprintf('Microstate %d', Labels(MSNum)));
        
        if isPeak
            timeInt = 1;
        else
            timeInt = 1 / FPS;
        end
        img{1, i} = frame2im(getframe(fig));
        [A,map] = rgb2ind(img{1, i},256);
        
        if i == 1
            imwrite(A,map,[DFPATH OUTPUTFILE],'gif','LoopCount',Inf,...
                'DelayTime', timeInt);
        else
            imwrite(A,map,[DFPATH OUTPUTFILE],'gif','WriteMode','append',...
                'DelayTime', timeInt);
        end
    end
    close;
    
    plot(1:90, allCor)
    hold on;
    for i = 10:10:90
        plot([i i], [-0.5 0.5], 'k--');
    end
    for i = 5:10:85
        text(i, 0.4, sprintf('%d', Labels((i + 5) / 10)));
    end
    legend({'ROI 1 Corr 1', 'ROI 2 Corr 2', 'ROI 3 Corr 3', 'ROI 4 Corr 4'});
    title(sprintf('Interval = %d ms, state Winlen = %d ms', INTERVAL, WINLEN(WINLENINDEX)));
end