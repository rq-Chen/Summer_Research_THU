%% topoCon.m - illustrate the significant connections between electrodes
%
% Ruiqi Chen, 2019.10.7

close all;

%% Parameters

DFPATH = '..\new\crq\Language\shantianfang\restRef\';
CLFILE = 'PLIClusteringData_Abs.mat';
CHFILE = '..\new\crq\rest\avrRef\msGFP.mat';
PROTHRES = 0.025;  % illustrate strongest connection
isDirected = 0;  % will calculate the abs val if set to 0
WINLEN = [10 20 45 100 200 700 2000];
WINLENINDEX = 2:7;
K_CLUSTERS = 3:5;
RADIUS = 0.5150;


%% Preparation

load([DFPATH CLFILE], 'ClusterData');

load(CHFILE, 'chanlocs');
allTheta = -cell2mat({chanlocs(:).theta}) + 90;
% originally clockwise and starting from y-axis
allR = cell2mat({chanlocs(:).radius});
allX = allR .* cosd(allTheta);
allY = allR .* sind(allTheta);
circleX = RADIUS .* cosd(0:5:360);
circleY = RADIUS .* sind(0:5:360);


%% Printing electrodes

figure;

for i = WINLENINDEX
    for j = K_CLUSTERS        
        Centers = ClusterData{1, i}.ClustData{1, j}.Centers;
        currWinlen = WINLEN(i);
        currClu = j;
        for l = 1:currClu
            if isDirected
                cofMat = squareform(Centers(l,:));
                tmp = sort(Centers(l,:));
            else
                cofMat = abs(squareform(Centers(l,:)));
                tmp = sort(abs(Centers(l,:)));
            end
            THRES = tmp(ceil((1 - PROTHRES) * size(tmp, 2)));
            hold on;
            for ll = 1:size(cofMat, 1)
                plot(circleX, circleY, 'k');
                plot(allX, allY, '.k', 'MarkerSize', 12);
            end
            for ll = 1:size(cofMat, 1)
                for lll = ll + 1 : size(cofMat, 1)
                    if cofMat(ll, lll) > THRES  % strong connection
                        tmp = cofMat(ll, lll);
                        if tmp < 0.1
                            colorTmp = [0 1 0];  % green
                        elseif tmp < 0.2
                            colorTmp = [0.8500 0.3250 0.0980];  % orange
                        else
                            colorTmp = [1 0 0];  % red
                        end
                        plot(allX([ll, lll]), allY([ll, lll]), ...
                            'Color', colorTmp, 'LineWidth', ...
                            2 * tmp);
                    end
                end
            end
            title(sprintf('State %d of %d', l, currClu));
            axis equal; axis off;
            saveas(gca, [DFPATH sprintf('PLIillustrate_Abs/%d/%d/%dTPCon.jpg', ...
                currWinlen, currClu, l)]);
            clf;
        end
    end
end
close all;