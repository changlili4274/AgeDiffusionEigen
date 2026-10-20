%**************************************************************************
    %> Filename: difmat_chebyshev_gauss_lobatto.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-01-28 15:04:30
    %> Description: Compute the Chebyshev-Gauss-Lobatto nodes, quadrature
    %  weights, and the first-order differentiation matrix on a given interval.
%**************************************************************************
function [x, w, DL] = difmat_chebyshev_gauss_lobatto(a, b, N)

%  INPUT:
%       a, b = extremes of the interval
%          M = number of nodes
%  OUTPUT:
%          x =  column vector containing N-Chebyshev extremal nodes in [a, b]
%          w =  column vector containing the weights for the Clenshaw-Curtis
%               quadrature associated to the nodes in x
%          DL = differentiation matrix relevant to x
%
% REFERENCES:
%            [1] Greg von Winckel (2021). Fast Clenshaw-Curtis Quadrature
%                (https://www.mathworks.com/matlabcentral/fileexchange/
%                6911-fast-clenshaw-curtis-quadrature), MATLAB Central
%                File Exchange. Retrieved November 16, 2021.
%
%            [2] L.N. Trefethen, "Spectral methods in Matlab", SIAM, 2000, 
%                DOI: 10.1137/1.9780898719598.
%
%            [3] J. A. Weideman and S. C. Reddy, "A MATLAB differentiation 
%                matrix suite", ACM Trans. Math. Software, 26 (2000), pp. 
%                465-519, https://doi.org/10.1145/365723.365727.

%%
if nargin < 3 || a >= b
    error('a must be less than b.');
end

if N < 2 || fix(N) ~= N
    error('N must be an integer >= 2 (number of nodes).');
end

bma = b - a;
if N == 2
    x = [a, b];
    w = 0.5*bma*[1, 1];
else
    n = N - 1;

    c = zeros(N, 2);
    c(1:2:N, 1) = (2./[1 1-(2:2:n).^2 ])';
    c(2, 2) = 1;

    f = real(ifft([c(1:N,:); c(n:-1:2,:)]));
    w = bma*([f(1,1); 2*f(2:n,1); f(N,1)])/2;
    x = 0.5*((b + a) + n*bma*f(1:N, 2));

    x = flip(x);
    w = flip(w);
end

c = [2; ones(N-2,1); 2].*(-1).^(0:N-1)';

X = repmat(x, 1, N);
dX = X - X';

DL = (c*(1./c)')./(dX + (eye(N)));
DL = DL - diag(sum(DL, 2));  

end
