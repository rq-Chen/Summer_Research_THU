%% illuSub.m - illustrate the functional connectivity change
%
% Ruiqi Chen, 2019.10.12
%
% Define "correspondent microstates" between different conditions as having
% the largest corrcoef. After that, the resting state (or its projection if
% specified) will be subtracted from the task state, and the remaining
% functional connectivity pattern will be illustrated.

clc; close all; clear;

%% Parameters

RESTFL = '../new/ydy/rest/avrRef/open/ClusteringData.mat';
TASKFL = '../new/ydy/Language/shantianfang/avrRef/ClusteringData.mat';
CHFILE = '../new/ydy/chanlocs.mat';
OUTPUTFD = '../new/ydy/cmpStates/';
OUTPUTFILE = [OUTPUTFD 'subClusteringData.mat'];
MINUSMODE = 0;  % 0:simply subtracting the resting state; 1:the projection
WINLEN = [10 20 45 100 200 700 2000];
WINLENINDEX = 2:7;
K_CLUSTERS = 3:5;
PROTHRES = 0.025;  % illustrate largest changes (top & bottom)
RADIUS = 0.5150;

%% Preparation

load(CHFILE, 'chanlocs');
allTheta = -cell2mat({chanlocs(:).theta}) + 90;
% originally clockwise and starting from y-axis
allR = cell2mat({chanlocs(:).radius});
allX = allR .* cosd(allTheta);
allY = allR .* sind(allTheta);
circleX = RADIUS .* cosd(0:5:360);
circleY = RADIUS .* sind(0:5:360);


%% Subtracting

load(RESTFL, 'ClusterData');
restClData = ClusterData;
load(TASKFL, 'ClusterData');
newClData = cell(length(WINLEN), 5);   % length() returns the size of the largest dim
figure;

for i = WINLENINDEX
    for j = K_CLUSTERS
        restCenter = restClData{1, i}.ClustData{1, j}.Centers;
        Center = ClusterData{1, i}.ClustData{1, j}.Centers;
        cofMat = corrcoef([restCenter' Center']);
        cofMat = cofMat(1:j, j + 1:end);  % rows: resting; columns: task
        % then choose a best permutation to get the largest sum of corrcoef
        allCorrSum = zeros(factorial(j), 1);
        allPer = perms(1:j);
        % each row of allPer represents a correspondence, e.g 1324 means
        % resting state - task state : 1-1 / 2-3 / 3-2 / 4-4
        for k = 1:factorial(j)
            for l = 1:j
                allCorrSum(k) = allCorrSum(k) + cofMat(l, allPer(k,l));
            end
        end
        [~, scInd] = max(allCorrSum);
        corres = allPer(scInd, :);
        
        SAVEFD = [OUTPUTFD sprintf('%d/%d/', WINLEN(i), j)];
        mkdir(SAVEFD);
        currData = zeros(j, length(Center));
        for k = 1:j
            
            if MINUSMODE == 0
                currData(k, :) = (Center(corres(k),:)) - (restCenter(k, :));
            else
                currData(k, :) = (Center(corres(k),:)) - (restCenter(k, :)) ...
                    * cosVal((Center(corres(k),:)), (restCenter(k, :)));
            end
            
            htMat = squareform(currData(k, :));
            heatmap(htMat, 'Colormap', parula);
            caxis([-0.5 0.5]); colorbar off;
            title(sprintf('task state %d - resting state %d, corr = %f', ...
                corres(k), k, cofMat(k, corres(k))));
            saveas(gca, [SAVEFD sprintf('%dsub%d.jpg', corres(k), k)]);
            clf;
            
            hold on;
            tmp = sort(currData(k,:), 'descend');
            THRES1 = tmp(ceil(PROTHRES * size(tmp, 2)));
            THRES2 = tmp(ceil((1 - PROTHRES) * size(tmp, 2)));
            for ll = 1:size(htMat, 1)
                plot(circleX, circleY, 'k');
                plot(allX, allY, '.k', 'MarkerSize', 12);
            end  
            for ll = 1:size(htMat, 1)
                for lll = ll + 1 : size(htMat, 1)
                    if htMat(ll, lll) > THRES1 || htMat(ll, lll) < THRES2
                        tmp = htMat(ll, lll);
                        % scale tmp to 1-64 for parula colormap
                        if abs(tmp) > 0.5
                            tmp = 0.5 * sign(tmp);
                        end
                        colormap = parula;
                        colorTmp = colormap(ceil((tmp + 0.5) * 63 + 1),:);
                        plot(allX([ll, lll]), allY([ll, lll]), ...
                            'Color', colorTmp, 'LineWidth', ...
                            2 * abs(htMat(ll, lll)));
                    end
                end
            end
            title(sprintf('task state %d - resting state %d, corr = %f', ...
                corres(k), k, cofMat(k, corres(k))));
            axis equal; axis off;
            colorbar('Ticks',[0, 0.2,0.4,0.6,0.8, 1],...
                'TickLabels',{'-0.5', '-0.3','-0.1','0.1','0.3', '0.5'})
            saveas(gca, [SAVEFD sprintf('%dsub%dTP.jpg', corres(k), k)]);
            clf;            
        end
        newClData{i, j} = currData;
    end
end

save(OUTPUTFILE, 'newClData');
close all;


function csVal = cosVal(A, B)
    lA = sqrt(dot(A, A));
    lB = sqrt(dot(B, B));
    dAB = dot(A, B);
    csVal = dAB / (lA * lB);
end









