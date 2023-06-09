function M = multmat2d(n, A, lambda_x, lambda_y)
%MULTMAT2D   2D multiplication matrix for the ultraspherical spectral method.
%   M = MULTMAT2D(N, A, LAMBDA_X, LAMBDA_Y) forms the N^2 x N^2
%   multiplication matrix in the tensor-product basis C^{(LAMBDA_Y)}(y) *
%   C^{(LAMBDA_x)}(x) for the function with bivariate Chebyshev
%   coefficients given by the matrix A.

if ( nargin < 3 ), lambda_x = 0;        end
if ( nargin < 4 ), lambda_y = lambda_x; end

if ( lambda_x == 0 && lambda_y == 0 )
    M = multmat2d_chebT(n, A);
elseif ( lambda_x >= 0 && lambda_y >= 0 )
    % Convert coefficients from Chebyshev to ultraspherical
    Sx = convertmat(n, 0, lambda_x-1);
    Sy = convertmat(n, 0, lambda_y-1);
    A = Sy * A * Sx.';
    if ( lambda_x > 0 && lambda_y > 0 )
        % Ultraspherical in x and y
        M = multmat2d_ultraS(n, A, lambda_x, lambda_y);
    elseif ( lambda_x == 0 )
        % Ultraspherical in y, Chebyshev in x
        M = multmat2d_ultraS_chebT(n, A, lambda_y);
    elseif ( lambda_y == 0 )
        % Chebyshev in y, ultraspherical in x
        M = multmat2d_chebT_ultraS(n, A, lambda_x);
    end
else
    error('Invalid ultraspherical parameter.')
end

end

function M = multmat2d_ultraS_chebT(n, A, lambda_y)
%MULTMAT2D_ULTRAS_CHEBT   Compute the 2D N^2 x N^2 multiplication matrix for the function
%
%   f(x,y) = sum_j sum_k  A(j,k) C^{(lambda_y)}_j(y) T_k(x)

m = 2*n;
Mold = speye(m, m);

Mx = spdiags(ones(m,2)/2,[-1 1], m, m);
Mx(2,1) = 1;

d1 = [1 2*lambda_y:2*lambda_y+m-2] ./ [1 2*((lambda_y+1):lambda_y+m-1)];
d2 = (1:m) ./ (2*(lambda_y:lambda_y+m-1));
B = [d2' zeros(m,1) d1'];
My = spdiags(B, [-1 0 1], m, m);
Mj = 2*lambda_y*My;

% j = 1:
M1 = multmat1d_chebT(n, A(1,:), Mx);
M = kron(M1, Mold(1:n,1:n));

for j = 2:n
    M1 = multmat1d_chebT(n, A(j,:), Mx);
    M = M + kron(M1, Mj(1:n,1:n));
    Mnew = (2*(j+lambda_y-1)/j)*My*Mj - ((j+2*lambda_y-2)/j)*Mold;
    Mold = Mj; Mj = Mnew;
end

end

function M = multmat2d_chebT_ultraS(n, A, lambda_x)
%MULTMAT2D_CHEBT_ULTRAS   Compute the 2D N^2 x N^2 multiplication matrix for the function
%
%   f(x,y) = sum_j sum_k  A(j,k) T_j(y)C^{(lambda_x)}_k(x)

m = 2*n;
Mold = speye(m, m);
My = spdiags(ones(m,2)/2,[-1 1], m, m);
My(2,1) = 1;
Mj = My;

d1 = [1 2*lambda_x:2*lambda_x+m-2] ./ [1 2*((lambda_x+1):lambda_x+m-1)];
d2 = (1:m) ./ (2*(lambda_x:lambda_x+m-1));
B = [d2' zeros(m,1) d1'];
Mx = spdiags(B, [-1 0 1], m, m);

% j = 1:
M1 = multmat1d_ultraS(n, A(1,:), Mx, lambda_x);
M = kron(M1, Mold(1:n,1:n));

for j = 2:n
    M1 = multmat1d_ultraS(n, A(j,:), Mx, lambda_x);
    M = M + kron(M1, Mj(1:n,1:n));
    Mnew = 2*My*Mj - Mold;
    Mold = Mj; Mj = Mnew;
end

end

function M = multmat2d_chebT(n, A)
%MULTMAT2D_CHEBT   Compute the 2D N^2 x N^2 multiplication matrix for the function
%
%   f(x,y) = sum_j sum_k  A(j,k) T_j(y)T_k(x)

m = 2*n;
Mold = speye(m, m);
Mx = spdiags(ones(m,2)/2,[-1 1], m, m);
Mx(2,1) = 1;
Mj = Mx;

% j = 1:
M1 = multmat1d_chebT(n, A(1,:), Mj);
M = kron(M1, Mold(1:n,1:n));

for j = 2:n
    M1 = multmat1d_chebT(n, A(j,:), Mx);
    M = M + kron(M1, Mj(1:n,1:n));
    Mnew = 2*Mx*Mj - Mold;
    Mold = Mj; Mj = Mnew;
end

end

function M = multmat2d_ultraS(n, A, lambda_x, lambda_y)
%MULTMAT2D_ULTRAS   Compute the 2D N^2 x N^2 multiplication matrix for the function
%
%   f(x,y) = sum_j sum_k  A(j,k) C^{(lambda_y)}_j(y)C^{(lambda_x)}_k(x)

m = 2*n;
Mold = speye(m, m);

d1 = [1 2*lambda_x:2*lambda_x+m-2] ./ [1 2*((lambda_x+1):lambda_x+m-1)];
d2 = (1:m) ./ (2*(lambda_x:lambda_x+m-1));
B = [d2' zeros(m,1) d1'];
Mx = spdiags(B, [-1 0 1], m, m);

d1 = [1 2*lambda_y:2*lambda_y+m-2] ./ [1 2*((lambda_y+1):lambda_y+m-1)];
d2 = (1:m) ./ (2*(lambda_y:lambda_y+m-1));
B = [d2' zeros(m,1) d1'];
My = spdiags(B, [-1 0 1], m, m);

Mj = 2*lambda_y*My;

% j = 1:
M1 = multmat1d_ultraS(n, A(1,:), Mx, lambda_x);
M = kron(M1, Mold(1:n,1:n));

for j = 2:n
    M1 = multmat1d_ultraS(n, A(j,:), Mx, lambda_x);
    M = M + kron(M1, Mj(1:n,1:n));
    Mnew = (2*(j+lambda_y-1)/j)*My*Mj - ((j+2*lambda_y-2)/j)*Mold;
    Mold = Mj; Mj = Mnew;
end

end
