function S = convertmat(n, k1, k2)
%CONVERTMAT   Conversion matrix for the ultraspherical spectral method.
%   S = CONVERTMAT(N, K1, K2) computes the N x N matrix realization of the
%   conversion operator between two bases of ultrapherical polynomials. The
%   matrix S maps N coefficients in a C^{(K1)} basis to N coefficients in a
%   C^{(K2 + 1)} basis, where C^{(K)} denotes ultraspherical polynomial
%   basis with parameter K. If K2 < K1, S is the N x N identity matrix.

S = speye(n);
for s = k1:k2
    S = spconvert(n, s) * S;
end

end

function T = spconvert(n, lam)
%SPCONVERT   Compute sparse representation for conversion operators.
%   SPCONVERT(N, LAM) returns the truncation of the operator that
%   transforms C^{(LAM)} (ultraspherical polynomials) to C^{(LAM+1)}. The
%   truncation gives back a matrix of size N x N.

% Relation is: C_n^(lam) = (lam/(n+lam))(C_n^(lam+1) - C_{n-2}^(lam+1))

if ( lam == 0 )
    dg = .5*ones(n-2,1);
    T = spdiags([1 0 ; .5 0 ; dg -dg], [0 2], n, n);
else
    dg = lam./(lam+(2:n-1))';
    T = spdiags([1 0 ; lam./(lam+1) 0 ; dg -dg], [0 2], n, n);
end

end
