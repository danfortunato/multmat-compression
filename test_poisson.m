rng(0)
u = randnfun2(1);
u = u - mean2(u);
f = lap(u);

dir = 1;
neu = 0;
[ux, uy] = grad(u);
g_w = dir*u(-1,:) - neu*ux(-1,:);
g_e = dir*u( 1,:) + neu*ux( 1,:);
g_s = dir*u(:,-1) - neu*uy(:,-1);
g_n = dir*u(:, 1) + neu*uy(:, 1);

n = 40;
method = 'ultraspherical';

v = poisson(n, f, g_w, g_e, g_s, g_n, dirichlet=dir, neumann=neu, method=method);
norm(u-v)
