%% cmpPLI.m - compute the signed Phase Lag Index over time
%
% Ruiqi Chen, 2019/10/25
%
% IMPORTANT: in order to get the original Phase Lag Index, use
%   abs(mean(cmpPLI(...))) to get the average value and ignore the
%   direction of coupling.
%
% Usage:
%   PL = cmpPLI(A, B); - two signal represented by column vectors
%       A, B and PL are column vectors of the same length
%   PL = cmpPLI(A); - n signal represented by each column of the matrix A
%       A is an m*n matrix and PL is an m * (n*(n-1)/2) matrix
%
% The Index:
%   Introduced in (Stam, 2007), doi: 10.1002/hbm.20346
%   Here we didn't do averaging or changing into absolute value. Therefore,
%     PL(t) begin positive indicates A is in advance relative to B at time
%     t, and vice versa.

function PL = cmpPLI(A, B)

if nargin == 0
    PL = [];
    return;
end

if nargin == 1
    if ~ismatrix(A)
        error('cmpPLI error: Dimension of the first input is not 2.');
    end
    if size(A, 1) ~= length(A)
        warnMsg = ['cmpPLI warning: A has more columns then rows! ' ...
            'Please check whether each column represents a sequence.'];
%         warning(warnMsg);
    end
    PL = nan(size(A, 1), nchoosek(size(A, 2), 2));
    loopI = 1;
    for i = 1:size(A, 2) - 1
        for j = i + 1:size(A, 2)
            PL(:, loopI) = cmpPLI(A(:, i), A(:, j));
            loopI = loopI + 1;
        end
    end
    return;
end

if nargin > 2
    error('cmpPLI error: too many inputs');
end

if ~iscolumn(A) || ~iscolumn(B) || length(A) ~= length(B)
    error(['cmpPLI error: either the input are not column vectors, ' ...
        'or they have different length']);
end

anaSig = hilbert([A B]);
allPhi = atan(imag(anaSig) ./ real(anaSig));
PL = sign(allPhi(:, 1) - allPhi(:, 2));

end