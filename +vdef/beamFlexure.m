function F = beamFlexure(x, h, E, opts)
%BEAMFLEXURE Tidal deflection of a floating elastic beam of varying thickness.
%   F = BEAMFLEXURE(x, h, E, opts) solves the Euler-Bernoulli beam equation
%   for ice in a grounding zone loaded by the hydrostatic pressure of the
%   ocean beneath it, after Holdsworth (1969):
%
%       d^2/dx^2 [ D(x) d^2w/dx^2 ] = rho_w * g * [ A0 - w(x) ]
%
%   with the flexural rigidity of the beam
%
%       D(x) = E h(x)^3 / (12 (1 - nu^2)).
%
%   w is the UPWARD deflection of the ice surface relative to its unstressed
%   position, A0 the far-field sea level, and the right-hand side the net
%   upward pressure that the difference between the two produces.
%
%   COORDINATES AND BOUNDARY CONDITIONS. x increases SEAWARD. The landward
%   boundary x(1) is the inward limit of tidal flexure, where the ice is
%   held against the bed:
%
%       w = 0  and  dw/dx = 0    at x(1)
%
%   The seaward boundary x(end) is out on the freely floating shelf, where
%   the ice has taken up the tide and gone flat:
%
%       w = A0  and  dw/dx = 0   at x(end)
%
%   This is Vaughan (1995)'s set. The seaward pair is a statement about the
%   far field, so x(end) must sit several flexural lengths seaward of the
%   hinge or it will artificially stiffen the beam; F.n_lambda reports how
%   many are actually there and is the number to check before believing a
%   result. The landward pair does NOT assume a rigid bed: any flexure of
%   the bed itself is absorbed into the effective modulus that a fit to
%   observed surface flexure returns, which is why that fitted quantity is
%   E* and not the Young's modulus of ice.
%
%   NUMERICS. Central differences on a uniform grid, after Jacquot and Dewey
%   (2001). The variable-rigidity operator is discretised in its self-adjoint
%   form,
%
%     (D w'')''|_i ~ [ D_{i-1} w''_{i-1} - 2 D_i w''_i + D_{i+1} w''_{i+1} ]/dx^2,
%
%   which gives the symmetric five-point stencil below and is second-order
%   accurate for smooth h. Rows 1, 2, N-1 and N carry the four boundary
%   conditions, the interior rows carry the equation.
%
%   Two things are done to keep the linear system solvable in double
%   precision. First the problem is nondimensionalised: lengths by the
%   reference flexural length lambda = (4 D_ref / (rho_w g))^(1/4) and
%   deflections by A0, which turns the equation into
%   (1/4) d^2/dxi^2[ d w'' ] + w = 1 with everything O(1). Second, each row
%   is divided by its own largest coefficient, so the boundary rows and the
%   interior rows are not left ~1e10 apart in scale. Without both, a
%   biharmonic operator on a few hundred nodes loses most of its digits.
%
%   INPUTS
%     x     N-vector of node positions [m], increasing, UNIFORMLY spaced
%     h     ice thickness [m]: scalar, or an N-vector on the same nodes
%     E     effective Young's modulus [Pa]
%     opts  all optional:
%       .nu     Poisson's ratio (default 0.3, as in the flexure literature)
%       .rho_w  seawater density [kg/m^3] (default vdef.constants)
%       .g      gravitational acceleration [m/s^2] (default vdef.constants)
%       .A0     far-field sea level [m] (default 1, giving a unit response
%               that scales linearly - the equation is linear in A0)
%
%   RETURNS
%     F.x, F.h    the grid and the thickness on it
%     F.w         N-vector of deflection [m], w(1) = 0 and w(end) = A0
%     F.dwdx      slope [-], second-order central differences
%     F.d2wdx2    curvature [1/m]. Bending strain is proportional to this:
%                 eps_xx(z) = -(z - z_n) * d2w/dx2 for a neutral surface at
%                 z_n, which is how a column-strain measurement rather than
%                 a surface-height measurement sees the same beam.
%     F.D         flexural rigidity [Pa m^3]
%     F.lambda    local flexural length (4D/(rho_w g))^(1/4) [m]
%     F.n_lambda  seaward domain length in flexural lengths - see above
%     F.residual  max |A*w - b| on the equilibrated system, a solve check
%
%   See also vdef.invertElasticModulus, scripts/diagnostics/elastic_modulus.m.

if nargin < 4 || isempty(opts), opts = struct(); end
C = vdef.constants();
if ~isfield(opts,'nu')    || isempty(opts.nu),    opts.nu    = 0.3;       end
if ~isfield(opts,'rho_w') || isempty(opts.rho_w), opts.rho_w = C.rho_sea; end
if ~isfield(opts,'g')     || isempty(opts.g),     opts.g     = C.g;       end
if ~isfield(opts,'A0')    || isempty(opts.A0),    opts.A0    = 1;         end

x = x(:);
N = numel(x);
assert(N >= 7, 'beamFlexure needs at least 7 nodes, got %d', N);

dx = diff(x);
assert(all(dx > 0), 'x must be strictly increasing (it runs landward to seaward)');
assert(max(abs(dx - dx(1))) <= 1e-6*dx(1), ...
  'x must be uniformly spaced; spacing varies by %.3g m', max(abs(dx - dx(1))));
dx = mean(dx);

if isscalar(h), h = repmat(h, N, 1); else, h = h(:); end
assert(numel(h) == N, 'h has %d values for %d nodes', numel(h), N);
assert(all(isfinite(h)) && all(h > 0), 'h must be finite and positive everywhere');
assert(isscalar(E) && isfinite(E) && E > 0, 'E must be a positive scalar [Pa]');

k = opts.rho_w * opts.g;                       % foundation modulus [Pa/m]
D = E * h.^3 / (12 * (1 - opts.nu^2));         % flexural rigidity [Pa m^3]

% Nondimensionalise. lambda is a single reference length for the whole
% beam, so the rigidity ratio d carries all of the thickness variation.
D_ref  = mean(D);
lambda = (4 * D_ref / k)^(1/4);
dxi    = dx / lambda;
d      = D / D_ref;

% Assemble (1/4) d^2/dxi^2[ d w'' ] + w = 1, both sides times 4*dxi^4 so the
% stencil coefficients come out O(d) instead of O(dxi^-4).
fnd = 4 * dxi^4;
b   = zeros(N, 1);

ii   = (3:N-2).';
rows = [ii; ii; ii; ii; ii];
cols = [ii-2; ii-1; ii; ii+1; ii+2];
vals = [ d(ii-1)
        -2*(d(ii-1) + d(ii))
         d(ii-1) + 4*d(ii) + d(ii+1) + fnd
        -2*(d(ii) + d(ii+1))
         d(ii+1) ];
b(ii) = fnd;

% Landward: w = 0, then dw/dx = 0 as a second-order one-sided difference.
% Seaward: the mirror of that slope condition, then w = 1.
rows = [rows; 1;   2;  2;  2;  N-1; N-1; N-1; N];
cols = [cols; 1;   1;  2;  3;  N-2; N-1; N;   N];
vals = [vals; 1;  -3;  4; -1;  1;  -4;   3;   1];
b([1 2 N-1 N]) = [0; 0; 0; 1];

A = sparse(rows, cols, vals, N, N);

% Row equilibration. The boundary rows are O(1) and the interior rows are
% O(d); scaling each row by its own largest entry leaves the solution
% unchanged and the condition number an order of magnitude better.
scale = full(max(abs(A), [], 2));
scale(scale == 0) = 1;
A = sparse(1:N, 1:N, 1./scale, N, N) * A;
b = b ./ scale;

W = A \ b;
assert(all(isfinite(W)), 'the beam solve did not converge to a finite deflection');

F = [];
F.x        = x;
F.h        = h;
F.w        = opts.A0 * W;
F.D        = D;
F.lambda   = (4 * D / k).^(1/4);
F.n_lambda = (x(end) - x(1)) / lambda;
F.residual = max(abs(A*W - b));

% Derivatives by central differences, one-sided at the two ends.
w = F.w;
F.dwdx    = zeros(N,1);
F.d2wdx2  = zeros(N,1);
F.dwdx(2:N-1)   = (w(3:N) - w(1:N-2)) / (2*dx);
F.dwdx(1)       = (-3*w(1) + 4*w(2) - w(3)) / (2*dx);
F.dwdx(N)       = ( 3*w(N) - 4*w(N-1) + w(N-2)) / (2*dx);
F.d2wdx2(2:N-1) = (w(3:N) - 2*w(2:N-1) + w(1:N-2)) / dx^2;
F.d2wdx2(1)     = F.d2wdx2(2);
F.d2wdx2(N)     = F.d2wdx2(N-1);

end
