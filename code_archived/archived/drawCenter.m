%% drawCenter - illustrate the Cluster Center

clear;clc;close;

%% Parameters

K_CLUSTERS = [3 4 5];
WINLENS = 2:7;  % index
DFPATH = '..\new\crq\rest\restref\';
FILENAME = 'PLIClusteringData.mat';
savePic = 1;


%% Generate Pictures

figure;
load([DFPATH FILENAME], 'ClusterData');
for i = WINLENS
    winlen = ClusterData{1, i}.Winlen;
    lag = ClusterData{1, i}.Lag;
    for j = K_CLUSTERS
        Pts = ClusterData{1, i}.ClustData{1, j}.Pt_Selected;
        Sil = ClusterData{1, i}.ClustData{1, j}.Sil_Selected;
        Labels = ClusterData{1, i}.ClustData{1, j}.Idx;
        Labels = Labels(Pts);
        fprintf('winlen = %d, lag = %d\n', winlen, lag);
        for l = 1:j
            sTitle = sprintf('winlen = %d, lag = %d, state %d of %d, Sil val = %f', ...
                winlen, lag, l, j, mean(Sil(Labels == l)));
            fprintf('%s\n', sTitle);
            
            heatMat = ClusterData{1, i}.ClustData{1, j}.Centers;
            heatMat = squareform(heatMat(l,:));            
%             for x = 1 : size(heatMat, 1)
%                 heatMat(x,x) = 1;
%             end
%             heatMat = tan(heatMat); % use tan to enlarge the difference
            heatmap(heatMat, 'Colormap', parula);
            caxis([-0.3 0.3]);
            title(sTitle);
            colorbar off;
            if savePic
                folderName = [DFPATH 'PLIillustrate/' num2str(winlen)...
                    '/' num2str(j) '/'];
                mkdir(folderName);
                saveas(gca, [folderName sprintf('%dlag%d_%d_%d.jpg',...
                    winlen, lag, l, j)]);
            end
            
        end
        
        %             MDSMat = ClusterData{1, i}.ClustData{1, j}.MDS_Selected;
        %             figure;
        %             gscatter(MDSMat(:,1), MDSMat(:,2), Labels);
        %             sTitle = sprintf('winlen = %d, lag = %d, k = %d, Sil val = %f', ...
        %                     winlen, lag, j, mean(Sil));
        %             title(sTitle);
        %             if savePic
        %                 folderName = [DFPATH 'illustrate/' num2str(winlen)...
        %                     '/' num2str(j) '/'];
        %                 mkdir(folderName);
        %                 saveas(gca, [folderName sprintf('%dlag%d.jpg',...
        %                     winlen, lag)]);
        %                 close all;
        %             end
        
    end
end
close all;
