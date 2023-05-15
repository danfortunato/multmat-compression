function [u, A, rhs] = poisson(n, f, g_w, g_e, g_s, g_n, opts)

arguments
    n
    f
    g_w = chebfun(0)
    g_e = chebfun(0)
    g_s = chebfun(0)
    g_n = chebfun(0)
    opts.dirichlet
    opts.neumann
    opts.method = 'ultraspherical'
end

isDirichlet = isfield(opts, 'dirichlet');
isNeumann   = isfield(opts, 'neumann');

if ( ~isDirichlet && ~isNeumann )
    % Default to Dirichlet boundary conditions
    opts.dirichlet = 1;
    opts.neumann   = 0;
elseif ( isDirichlet && ~isNeumann )
    opts.neumann = 0;
elseif ( ~isDirichlet && isNeumann )
    opts.dirichlet = 0;
end

if ( opts.dirichlet == 0 && opts.neumann == 0 )
    error('No boundary condition specified.');
end

a = opts.dirichlet;
b = opts.neumann;

switch lower(opts.method)
    case 'ultraspherical'
        basis = chebpoly(0:n-1);
        f = coeffs2(f, n, n);
        g_w = chebcoeffs(g_w, n);
        g_e = chebcoeffs(g_e, n);
        g_s = chebcoeffs(g_s, n);
        g_n = chebcoeffs(g_n, n);
        S02 = ultraS.convertmat(n, 0, 1);
        D2  = ultraS.diffmat(n, 2);
        ipde = false(n);    ipde(1:n-2,1:n-2) = true;
        ibc  = false(n, 1); ibc(1:n-2) = true;
        chebargs = {'coeffs'};
    case 'collocation'
        basis = chebfun.lagrange(chebpts(n));
        f = sample(f, n, n);
        g_w = sample(g_w, n);
        g_e = sample(g_e, n);
        g_s = sample(g_s, n);
        g_n = sample(g_n, n);
        S02 = speye(n);
        D2 = diffmat(n, 2);
        ipde = false(n);    ipde(2:n-1,2:n-1) = true;
        ibc  = false(n, 1); ibc(2:n-1) = true;
        chebargs = {};
end

% Differential operator
L = kron(D2, S02) + kron(S02, D2);
f = S02 * f * S02.';

% Boundary conditions
I = speye(n);
dbasis = diff(basis);
basis_m = a*basis(-1) - b*dbasis(-1);
basis_p = a*basis( 1) + b*dbasis( 1);
bc_w = kron(basis_m, I);
bc_e = kron(basis_p, I);
bc_s = kron(I, basis_m);
bc_n = kron(I, basis_p);

% Corner compatibility conditions
west  = @(x) basis_m*x;
east  = @(x) basis_p*x;
south = @(y) basis_m*y;
north = @(y) basis_p*y;
corner_sw = south(bc_w) + west(bc_s);
corner_se = south(bc_e) + east(bc_s);
corner_nw = north(bc_w) + west(bc_n);
corner_ne = north(bc_e) + east(bc_n);

% Assemble the linear system
A = [ L(ipde,:)   ;  % PDE
      bc_w(ibc,:) ;  % Boundary conditions
      bc_e(ibc,:) ;  % .
      bc_s(ibc,:) ;  % .
      bc_n(ibc,:) ;  % .
      corner_sw   ;  % Corner compatibility conditions
      corner_se   ;  % .
      corner_nw   ;  % .
      corner_ne   ]; % .

rhs = [ f(ipde)                ;
        g_w(ibc)               ;
        g_e(ibc)               ;
        g_s(ibc)               ;
        g_n(ibc)               ;
        south(g_w) + west(g_s) ;
        south(g_e) + east(g_s) ;
        north(g_w) + west(g_n) ;
        north(g_e) + east(g_n) ];

% Add the integral constraint for the pure Neumann problem
if ( a == 0 )
    w = sum(basis);
    ww = kron(w, w);
    A = [ A ww.' ; ww 0 ];
    rhs = [ rhs ; 0 ];
end

% Solve
uu = A \ rhs;
U = reshape(uu(1:n^2), n, n);
u = chebfun2(U, chebargs{:});

end
