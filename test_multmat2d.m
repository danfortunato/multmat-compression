f = @(x,y) sin(x.*y);
g = @(x,y) exp(y).*x.^2;
h = @(x,y) f(x,y).*g(x,y);
n = 30;
tol = eps;

% Compute the bivariate Chebyshev coefficients
[xx, yy] = chebpts2(n);
cf = chebvals2chebcoeffs( chebvals2chebcoeffs(f(xx,yy)).' ).';
cg = chebvals2chebcoeffs( chebvals2chebcoeffs(g(xx,yy)).' ).';
ch = chebvals2chebcoeffs( chebvals2chebcoeffs(h(xx,yy)).' ).';
cf(abs(cf) < tol) = 0;
cg(abs(cg) < tol) = 0;

%% Test Chebyshev multiplication
Mf = multmat2d(n, cf, 0);
Mg = multmat2d(n, cg, 0);

cfg = Mf * cg(:); cfg = reshape(cfg, n, n);
cgf = Mg * cf(:); cgf = reshape(cgf, n, n);

norm(cfg - ch)
norm(cgf - ch)

%% Test ultraspherical C^{(1)}(y) C^{(2)}(x) multiplication
lambda_x = 2;
lambda_y = 1;
Mf = multmat2d(n, cf, lambda_x, lambda_y);
Mg = multmat2d(n, cg, lambda_x, lambda_y);

Sx = convertmat(n, 0, lambda_x-1);
Sy = convertmat(n, 0, lambda_y-1);
cf_lam = Sy * cf * Sx.';
cg_lam = Sy * cg * Sx.';

cfg_lam = Mf * cg_lam(:); cfg_lam = reshape(cfg_lam, n, n);
cgf_lam = Mg * cf_lam(:); cgf_lam = reshape(cgf_lam, n, n);

norm(cfg_lam - Sy * ch * Sx.')
norm(cgf_lam - Sy * ch * Sx.')
