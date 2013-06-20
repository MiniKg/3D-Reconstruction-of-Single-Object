function [S, good, ind] = screen_hao(x, y)

sze = size(x);
nframes = sze(2);
nfeatures = sze(1);

%compute feature location centroid of each frame and 
%recenter the 2D coordinates of feature points
phati = zeros(nframes, 2);
xrecentered = x;
yrecentered = y;
for j = 1 : nframes
    phati(j, 1) = sum(x(:, j)) / nfeatures;
    phati(j, 2) = sum(y(:, j)) / nfeatures;
    xrecentered(:, j) = x(:, j) - phati(j, 1);
    yrecentered(:, j) = y(:, j) - phati(j, 2);
end

%construct the 2m*n measurement matrix W
tempW = zeros(nfeatures, nframes*2);
tempW(:, 1:2:nframes*2) = xrecentered;
tempW(:, 2:2:nframes*2) = yrecentered;
W = tempW';

%perform the factorization into shape and motion
[U, S, V] = svd(W);
Uhat = U(:, 1:3);
Vhat = V(:, 1:3);
Shat = S(1:3, 1:3);
Msvd = Uhat * sqrt(Shat);
Ssvd = sqrt(Shat) * Vhat';
newtracks = Msvd * Ssvd;

%discard tracks where the maximum deviation from the measured locations is
%too large
deviation = W - newtracks;
[~, bad2y] = find(deviation > 1.0);
bad2ynew = unique(bad2y);
W_size = size(W);
good = (1:W_size(2));
good = good(:);
good(bad2ynew, :) = [];
size(W)
size(bad2ynew)
W(:, bad2ynew) = [];   %discard corresponding columns

%rearrange the measurement matrix for function orthonormalize
W_sze = size(W);
newW = W;
newW(1:nframes, :) = W(1:2:W_sze(1), :);
newW(nframes+1:W_sze(1), :) = W(2:2:W_sze(1), :);

[U, S, V] = svd(newW);
Uhat = U(:, 1:3);
Vhat = V(:, 1:3);
Shat = S(1:3, 1:3);
Msvd = Uhat * sqrt(Shat);
Ssvd = sqrt(Shat) * Vhat';

[M, S] = orthonormalize(Msvd, Ssvd);
ind = good';
good = x(good', :);
end

