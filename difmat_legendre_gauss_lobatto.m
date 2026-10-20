%**************************************************************************
    %> Filename: difmat_legendre_gauss_lobatto.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-04-19 10:12:33
    %> Description: Compute the Legendre-Gauss-Lobatto nodes, quadrature
    %  weights, and the first-order differentiation matrix on a given interval.
%**************************************************************************
function [x, w, D] = difmat_legendre_gauss_lobatto(a, b, N)
%  INPUT:
%       a, b = extremes of the interval
%          N = number of nodes
%  OUTPUT:
%          x =  column vector containing N Legendre-Guass-Lobatto nodes 
%               on [a, b]
%          w =  column vector containing the weights for the quadrature 
%               associated to the nodes in x
%          D = differentiation matrix relevant to x

%%
if nargin < 3 || a >= b
    error('a must be less than b.');
end

if N < 2 || fix(N) ~= N
    error('N must be an integer >= 2 (number of nodes).');
end

% Compute LGL nodes, weights and Vandermonde matrix on [-1,1] (ascending order)
[x_ref, w_ref, P] = lglnodes(N);

% Construct first-order differentiation matrix on [-1,1]
D_ref = zeros(N);
P_last = P(:, N);   % values of P_{N-1}(x_i) at each node

for i = 1:N
    for j = 1:N
        if i ~= j
            D_ref(i, j) = (P_last(i) / P_last(j)) / (x_ref(i) - x_ref(j));
        end
    end
    D_ref(i, i) = -sum(D_ref(i, :));
end

%% Affine transformation from [-1,1] to [a,b]
x = (b - a)/2 * x_ref + (a + b)/2;
w = (b - a)/2 * w_ref;
D = (2 / (b - a)) * D_ref;

end


%% compute LGL nodes, weights and Legendre Vandermonde matrix on [-1,1]
function [x, w, P] = lglnodes(N)
%  INPUT:
%          N = number of nodes
%  OUTPUT:
%          x =  column vector containing N Legendre-Guass-Lobatto nodes 
%               on [-1, 1]
%          w =  column vector containing the weights for the quadrature 
%               associated to the nodes in x
%          P = Legendre Vandermonde matrix
%
% REFERENCES:
%            [1] C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, 
%                "Spectral Methods in Fluid Dynamics," Section 2.3. 
%                Springer-Verlag 1987.
%
%            [2] Greg von Winckel - 04/17/2004 Contact: gregvw@chtm.unm.edu 
%                (https://ww2.mathworks.cn/matlabcentral/fileexchange/4775-
%                legende-gauss-lobatto-nodes-and-weights?s_tid=srchtitle_
%                support_results_3_Gauss+Lobatto)
%
%            [3] J. A. Weideman and S. C. Reddy, "A MATLAB differentiation 
%                matrix suite", ACM Trans. Math. Software, 26 (2000), pp. 
%                465-519, https://doi.org/10.1145/365723.365727.

M = N - 1; % polynomial degree

% Initial guess: Chebyshev-Gauss-Lobatto nodes (descending order)
x = cos(pi*(0:M)'/M);

x_old = 2;
while max(abs(x - x_old)) > eps
    x_old = x;

    % Build Legendre Vandermonde matrix via recurrence
    P = zeros(N, N);
    P(:, 1) = 1;
    if N > 1
        P(:, 2) = x;
    end

    for k = 2:M
        P(:, k+1) = ((2*k-1)*x.*P(:, k) - (k-1)*P(:, k-1))/k;
    end

    % Newton update
    x = x_old - (x_old.*P(:, N) - P(:, M))./(N*P(:, N));
end

% Compute weights (using final P matrix)
w = 2./(M*N*(P(:, N).^2));

% Convert from descending order to ascending order
x = flip(x);
w = flip(w);
P = flip(P, 1);   % reverse rows (keep columns unchanged)

end
