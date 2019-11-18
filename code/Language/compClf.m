%% compClf - prepare for classifier training, executed after corrLanClust.m

ptEachWin = fix(SAMPLERATE * Winlen / 1000);  % number of timepoints within each window
ptStart = 1 + Lag;  % discard the data before it
clNum = fix((size(eegdata, 2) - ptStart + 1) ...
    / ptEachWin);  % number of points for clustering

svmInput = zeros(clNum, nchoosek(size(eegdata,1), 2) + 1);
for i = 1:clNum
    currentData = eegdata(:, (Lag + 1 + (i - 1) * Winlen) : ...
        (Lag + i * Winlen));  % get the data within this window
    cofMat = corrcoef(currentData');
    svmInput(i, 1:end - 1) = squareform(tril(cofMat,-1));
    % squareform doesn't check for symmetry, but does check for zero diag
end

LanIdx = zeros(1, size(eegdata, 2));
for i = 1:size(LanClustLabel, 1)
    startIdx = ceil(LanClustLabel{i, 1} * SAMPLERATE);
    endIdx = ceil(LanClustLabel{i, 2} * SAMPLERATE);
    if endIdx > size(eegdata, 2)
        LanClustLabel = LanClustLabel(1:i - 1, :);
        break;
    end
%     if LanClustLabel{i, 4} == ""
%         tmp = 3;
%     else
%         if LanClustLabel{i, 4} == "v"
%             tmp = 1;
%         else
%             tmp = 2;
%         end
%     end
%     LanIdx(startIdx:endIdx) = ones(1, endIdx - startIdx + 1) * tmp;

    LanIdx(startIdx:endIdx) = ones(1, endIdx - startIdx + 1);

end

Labels = zeros(1, clNum);
for i = 1:clNum
    svmInput(i, end) = max(LanIdx(Lag + 1 + (i - 1) * Winlen : Lag + i * Winlen));
end