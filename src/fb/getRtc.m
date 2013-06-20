%get the transformation parameters from two point clouds X and Y
%rotation R, translation t, scaling c
function [R, t, c] = getRtc(X, Y)

%get related parameters
cloud_size = size(X);
ux = sum(X, 2) / cloud_size(2);
uy = sum(Y, 2) / cloud_size(2);
Xdiff = X-repmat(ux, [1, cloud_size(2)]);
sigma2x = sum(sum(Xdiff.^2)) / cloud_size(2);
Ydiff = Y-repmat(uy, [1, cloud_size(2)]);
sigma2y = sum(sum(Ydiff.^2)) / cloud_size(2);
capsigmaxy = Ydiff * Xdiff' / cloud_size(2);

[U, D, V] = svd(capsigmaxy);
%only if det(capsigmaxy) > 0
S = eye(size(capsigmaxy));

%compute R, t, c
R = U * S * V';
traceDS = sum(diag(D*S));
c = 1/sigma2x * traceDS;
t = uy - c*R*ux;

end