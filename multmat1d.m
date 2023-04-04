function M = multmat1d(n, a, lambda)
%MULTMAT1D   1D multiplication matrix for the ultraspherical spectral method.
%   M = MULTMAT1D(N, A, LAMBDA) forms the N x N multiplication matrix in
%   the C^{(LAMBDA)} basis for the function with Chebyshev coefficients
%   given by the vector A.

if ( nargin < 3 ), lambda = 0; end

if ( isempty(a) )
    M = sparse(n, n);
    return
elseif ( numel(a) == 1 )
    M = a * speye(n);
    return
end

% Prolong or truncate coefficients
a = a(:);
if ( numel(a) < n )
    a = [a ; zeros(n-numel(a), 1)];
else
    a = a(1:n);
end

% Convert coefficients to the right basis
a = convertmat(n, 0, lambda-1) * a;

if ( isempty(a) )
    M = sparse(n, n);
else
    m = 2*n;
    if ( lambda == 0 )
        Mx = spdiags(ones(m,2)/2, [-1 1], m, m);
        Mx(2,1) = 1;
        M = multmat1d_chebT(n, a, Mx);
    elseif ( lambda > 0 )
        d1 = [1 2*lambda:2*lambda+m-2] ./ [1 2*((lambda+1):lambda+m-1)];
        d2 = (1:m) ./ (2*(lambda:lambda+m-1));
        B = [d2' zeros(m,1) d1'];
        Mx = spdiags(B, [-1 0 1], m, m);
        M = multmat1d_ultraS(n, a, Mx, lambda);
    else
        error('Invalid ultraspherical parameter.')
    end
end

end
