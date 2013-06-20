function [R3, S3] = orthonormalize(Rhat, Shat)
% [R3, S3] = orthonormalize(Rhat, Shat)
%
% Given the 2Fx3 matrix Rhat and the 3xP matrix Shat,
% this function returns an equivalent decomposition
% R3, S3 such that R3 is orthonormal.


[F,D] = size(Rhat); F = F/2;

% Build matrix Q such that Q * v = [1,...,1,0,...,0] where v is a six
% element vector containg all six distinct elements of the Matrix C

clear Q
for f = 1:F,
  g = f + F;
  h = g + F;
  Q(f,:) = zt2(Rhat(f,:), Rhat(f,:));
  Q(g,:) = zt2(Rhat(g,:), Rhat(g,:));
  Q(h,:) = zt2(Rhat(f,:), Rhat(g,:));
end

% Solve for v
rhs = [ones(2*F,1); zeros(F,1)];
v = Q \ rhs;

% C is a symmetric 3x3 matrix such that C = G * transpose(G)
n = 1;
for x = 1:D,
  for y = x:D,
    C(x,y) = v(n); C(y,x) = v(n);
    n = n+1;
  end;
end;

e = eig(C);
% disp(e);

if (any(e<= 0)),
  [uu,dd,vv] = svd(C);
  G = uu*sqrt(dd);
  cc = G*G';
  err = sum((cc(:)-C(:)).^2);
  if err>0.03,
    G = [];
    dbstack; keyboard
  end;
else
  G = sqrtm(C);
end

R3 = Rhat*G;
S3 = inv(G)*Shat;


%-------------------------------------------------

function M = zt2(i,j)

  D = length(i);
  M = [];
  for x = 1:D,
    for y = x:D,
      if x==y, M = [M, i(x)*j(y)]; 
      else,    M = [M, i(x)*j(y)+i(y)*j(x)];
      end;
    end;
  end;

