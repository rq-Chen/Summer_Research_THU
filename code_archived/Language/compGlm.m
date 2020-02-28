%% compGlm - prepare for glm classifing, executed after corrLanClust

N_TIMELEN = 31;
N_ST = -200;

glmInput = zeros(size(Corr, 2), 1 + 3 * N_TIMELEN);
glmInput(:, end) = Corr(CL_OI, :);  % correaltion with the state-of-interest

MUSICFILE = '13.wav';
[audioMat, fs] = audioread([DFPATH MUSICFILE]);
iCenter = 1 - N_ST / Winlen;

for i = 1:size(glmInput, 1)
    startSp = 1 + fix((i - 1) * Winlen / 1000 * fs);
    endSp = fix(i * Winlen / 1000 * fs);
    tmp = audioMat(startSp:endSp, 1) .^ 2;
    glmInput(i, iCenter) = mean(tmp);  % volume
end

for i = 1:size(LanClustLabel, 1)
    glmSIdx = 1 + fix(MusicLabel{i, 1} * 1000 / Winlen);
    glmEIdx = 1 + fix(MusicLabel{i, 2} * 1000 / Winlen);
    glmInput(glmSIdx, iCenter + 2 * N_TIMELEN) = MusicLabel{i, end};  % acoustic onset
    glmInput(glmSIdx, iCenter + N_TIMELEN) = 1;  % word onset
end

for i = 1:N_TIMELEN
    
    if i == iCenter
        continue;
    end
    startSp = i - iCenter;
    if startSp < 0
        glmInput(1: end + startSp, i:N_TIMELEN:i + 2 * N_TIMELEN) = ...
            glmInput(1 - startSp:end, iCenter:N_TIMELEN:iCenter + 2 * N_TIMELEN);
    else
        glmInput(1 + startSp:end, i:N_TIMELEN:i + 2 * N_TIMELEN) = ...
            glmInput(1: end - startSp, iCenter:N_TIMELEN:iCenter + 2 * N_TIMELEN);
    end        
    
end