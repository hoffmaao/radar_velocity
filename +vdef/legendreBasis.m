function [P, dP] = legendreBasis(x, K)
%LEGENDREBASIS Legendre polynomials and derivatives up to order K.
%   [P, dP] = LEGENDREBASIS(x, K) returns N x (K+1) matrices whose columns
%   are P_0..P_K and dP_0/dx..dP_K/dx evaluated at the N points x, using
%   the standard recurrences
%
%     (k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}
%     P'_{k+1}      = (2k+1) P_k + P'_{k-1}
%
%   x is expected on [-1, 1]. Legendre polynomials are used rather than a
%   raw monomial basis because they are orthogonal on [-1, 1], so with
%   roughly uniform depth sampling the normal matrix is close to diagonal
%   and the low-order coefficients (which carry the physics) are not
%   traded off against the high-order ones by the fit.
%
%   See also vdef.invertStrainRate.

x = x(:);
N = numel(x);
P  = zeros(N, K+1);
dP = zeros(N, K+1);

P(:,1)  = 1;
dP(:,1) = 0;
if K >= 1
  P(:,2)  = x;
  dP(:,2) = 1;
end
for k = 1:K-1
  P(:,k+2)  = ((2*k+1) * x .* P(:,k+1) - k * P(:,k)) / (k+1);
  dP(:,k+2) = (2*k+1) * P(:,k+1) + dP(:,k);
end

end
