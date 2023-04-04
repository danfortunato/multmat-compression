f = @(x) sin(3*x);
g = @(x) x.^2.*exp(x);
h = @(x) f(x).*g(x);
n = 30;
tol = eps;

% Compute the Chebyshev coefficients
xx = chebpts(n);
cf = chebvals2chebcoeffs(f(xx));
cg = chebvals2chebcoeffs(g(xx));
ch = chebvals2chebcoeffs(h(xx));
cf(abs(cf) < tol) = 0;
cg(abs(cg) < tol) = 0;

%% Test Chebyshev multiplication
Mf = multmat1d(n, cf, 0);
Mg = multmat1d(n, cg, 0);

cfg = Mf * cg;
cgf = Mg * cf;

norm(cfg - ch)
norm(cgf - ch)

%% Test ultraspherical C^{(2)} multiplication
lambda = 2;
Mf = multmat1d(n, cf, lambda);
Mg = multmat1d(n, cg, lambda);

S = convertmat(n, 0, lambda-1);
cf_lam = S * cf;
cg_lam = S * cg;

cfg_lam = Mf * cg_lam;
cgf_lam = Mg * cf_lam;

norm(cfg_lam - S * ch)
norm(cgf_lam - S * ch)
