function [xx, yy] = chebpts2(nx, ny, dom)
%CHEBPTS2   Chebyshev tensor product grid.
%   [XX, YY] = CHEBPTS2(N) constructs an N x N grid of tensor-product
%   Chebyshev points of the 2nd kind on [-1,1]^2.
%
%   [XX, YY] = CHEBPTS2(NX, NY) constructs an NX x NY grid of
%   tensor-product Chebyshev points of the 2nd kind on [-1,1]^2.
%
%   [XX, YY] = CHEBPTS2(NX, NY, DOM) constructs an NX x NY grid of
%   tensor-product Chebyshev points of the 2nd kind on the rectangle
%   [a,b] x [c,d], where DOM = [a b c d].
%
%   See also CHEBPTS.

if ( nargin < 2 ), ny = nx;           end
if ( nargin < 3 ), dom = [-1 1 -1 1]; end

x = chebpts(nx, dom(1:2));
y = chebpts(ny, dom(3:4));
[xx, yy] = meshgrid(x, y);

end
