%**************************************************************************
    %> Filename: lagrange_basis_barycentric.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-04-16 09:59:45
    %> Description: Values of the Lagrange basis functions in barycentric
    %  form at a given point.
%**************************************************************************
function basis_val = lagrange_basis_barycentric(x, xk, wk)
% Input arguments:
%   x: interpolation point
%   xk: interpolation nodes
%   wk: barycentric weights
%
% Output argument:
%   basis_val: values of all Lagrange basis functions at the point x (column vector)

N = length(xk);
basis_val = zeros(N, 1);

[~, idx] = min(abs(x - xk));
if abs(x - xk(idx)) < eps
    basis_val(idx) = 1;
    return;
end

basis_val = wk ./ (x - xk);
basis_val = basis_val / sum(basis_val);

end
