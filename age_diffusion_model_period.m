%**************************************************************************
    %> Filename: age_diffusion_model_period.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-01-29 19:10:33
    %> Description: Compute the eigenvalues and eigenfunctions of the
    %  age-structured local diffusion model with spatially periodic boundaries.
%**************************************************************************
function [eigenvalues, eigenfunctions] = age_diffusion_model_period(...
    a, s, r, T, d, beta, pi0, mu1, N, M, p, discontinuity_age, foldername)
% Input arguments:
%   a: age interval [a(1), a(end)]
%   s: left endpoint of the spatial interval; the spatial interval is [s, s+1]
%   lam: parameter of the weighted eigenvalue problem
%   T: spatial period
%   d: diffusion coefficient
%   beta: birth-rate function beta = beta(a, x)
%   pi0: natural survival probability pi0 = exp(-int_{0}^{a} mu0(s) ds), 
%        where mu0 is the natural mortality function
%   mu1: mortality function mu1 = mu1(a, x) depending on the spatial location
%   N: number of Chebyshev-Gauss-Lobatto nodes in the age direction
%   M: number of nodes in the spatial direction
%   p: plotting flag (plot if nonzero)
%   discontinuity_age: optional vector of discontinuities of the birth 
%      rate beta; if provided, piecewise adaptive integration is used
%   foldername: folder name, including the path, in which the data are saved
%
% Output arguments:
%   eigenvalues: approximate eigenvalues
%   eigenfunctions: approximate eigenfunctions

%%
% Obtain the Chebyshev-Gauss-Lobatto nodes and the Chebyshev 
% differentiation matrix in the age direction
[a_points, ~, a_DL] = difmat_chebyshev_gauss_lobatto(a(1), a(end), N);
% Obtain the nodes and the Fourier differentiation matrix 
% in the spatial direction
[s_points, D1, D2, ~] = difmat_fourier(s, s+1, M);

% Values of the known functions at the nodes
pi0_a_points = pi0(a_points);

[s_mesh, a_mesh] = meshgrid(s_points, a_points);
mu1_as_points = mu1(a_mesh, s_mesh);

%% Compute the integral coefficients of the boundary condition
% Barycentric weights of the Chebyshev-Gauss-Lobatto interpolation nodes
bary_weights = (-1).^(0:N-1)';
bary_weights(1) = 0.5*bary_weights(1);
bary_weights(end) = 0.5*bary_weights(end);

% Lagrange basis functions l_k(a)
lagrange_basis = @(a) lagrange_basis_barycentric(a, a_points, bary_weights);

if nargin < 12
    % The function of age has no discontinuities
    discontinuity_age = [];
end

% Integration settings
opts = {'RelTol', eps, 'AbsTol', eps, 'ArrayValued', true, ...
    'Waypoints', discontinuity_age};
integral_bc = zeros(N, M);

fprintf('Computing the age-integral boundary condition...');
for j = 1:M
    xj = s_points(j);
    f = @(a) beta(a, xj).*pi0(a).*lagrange_basis(a);
    integral_bc(:, j) = integral(f, a(1), a(end), opts{:});
end
fprintf('done.\n');

%% Assemble the global matrix
D = kron(speye(M), a_DL) ... 
    - kron(d/T^2*D2, speye(N)) ...
    + kron(2*d*r/T*D1, speye(N)) ... 
    + spdiags(mu1_as_points(:) - d*r^2, 0, N*M, N*M);

%% Modify the global matrix using the boundary condition in the age direction
j = 1;
for i = 1:N:N*M
    D(i, :) = zeros(1, N*M);
    D(i, i:i+N-1) = - integral_bc(:, j)';
    D(i, i) = D(i, i) + 1;
    j = j + 1;
end

B = kron(speye(M), spdiags([0; ones(N-1, 1)], 0, N, N));

%% Compute the eigenvalues and eigenfunctions
[eigenfunctions, eigenvalues] = eig(full(D\B));
eigenvalues = 1./diag(eigenvalues);

% Remove invalid eigenvalues and the corresponding eigenvectors
noninf_ind = find(abs(eigenvalues)~=Inf);
eigenvalues = eigenvalues(noninf_ind);
eigenfunctions = eigenfunctions(:, noninf_ind);

[~, ind] = sort(eigenvalues, 'ComparisonMethod', 'real');
eigenvalues = - eigenvalues(ind);
eigenfunctions = eigenfunctions(:, ind);

eigenfunctions = kron(ones(M, 1), pi0_a_points).*eigenfunctions;

%% Save data related to the principal eigenvalues and eigenfunctions
lambda1 = eigenvalues(1);
lambda3 = eigenvalues(3);
lambda5 = eigenvalues(5);

phi1 = eigenfunctions(:, 1);
phi3 = eigenfunctions(:, 3);
phi5 = eigenfunctions(:, 5);

% Print the principal eigenvalue
fprintf('λ1 = %.10f + %.10f i\n', real(lambda1), imag(lambda1));

% Folder name for saving data
if nargin < 13 || isempty(foldername)
    foldername = '.\age_diffusion_period_data\';
end

% Create the folder for saving data
if ~isfolder(foldername) 
    mkdir(foldername);
end

% Compute c = H(r, T)/r
c = lambda1/r;

filename = [foldername, 'N', num2str(N), '_M', num2str(M), ...
    '_r', num2str(r), '_T', num2str(T), '_d', num2str(d), '.mat'];
save(filename, 'r', 'T', 'd', 'c', 'a_points', 's_points', ...
    'lambda1', 'phi1', 'lambda3', 'phi3', 'lambda5', 'phi5');

%% Visualization
if p
    figure(1);

    % Plot the eigenvalues
    plot(real(eigenvalues), imag(eigenvalues), 'ko', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 6);
    xlabel('Re(\lambda)');
    ylabel('Im(\lambda)');
    % title('Eigenvalues');
    axis([-70, lambda1+3.0, -100, 100]);
    grid on;
    hold on;

    % Mark the principal eigenvalues
    plot(real(lambda1), imag(lambda1), 'ro', 'MarkerSize', 10, ...
        'LineWidth', 2, 'MarkerFaceColor', 'r');
    hold on;
    plot(real(lambda3), imag(lambda3), 'gs', 'MarkerSize', 10, ...
        'LineWidth', 2, 'MarkerFaceColor', 'g');
    hold on;
    plot(real(lambda5), imag(lambda5), 'bd', 'MarkerSize', 10, ...
        'LineWidth', 2, 'MarkerFaceColor', 'b');
    hold off;

    legend('', '\lambda_1', '\lambda_3', '\lambda_5', 'Location', 'west');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
    set(get(gca,'Children'), 'linewidth', 1.0);
    set(get(gca,'YLabel'), 'FontSize', 25, 'FontName', 'Times New Roman');
    set(get(gca,'XLabel'), 'FontSize', 25, 'FontName', 'Times New Roman');

    % Plot the principal eigenfunction
    figure(2);
    
    phi1 = reshape(phi1, N, M);
    phi1 = [phi1, phi1(:, 1)];
    surf([s_mesh, s_mesh(:, 1)+1], [a_mesh, a_mesh(:, 1)], abs(phi1));
    
    xlabel('\it{x}');
    ylabel('\it{a}');
    zlabel('\Phi_1');
    shading interp;
    grid on;
    % title('The principal eigenfunction');

    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
    set(get(gca,'Children'), 'linewidth', 1.0);
    set(get(gca,'YLabel'), 'FontSize', 25, 'FontName', 'Times New Roman');
    set(get(gca,'XLabel'), 'FontSize', 25, 'FontName', 'Times New Roman');
end

end
