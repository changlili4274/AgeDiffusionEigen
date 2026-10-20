%**************************************************************************
    %> Filename: example1.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-02-02 17:33:54
    %> Description: Numerical example of the eigenvalue problem for the
    %  age-structured local diffusion model with nonperiodic boundaries.
%**************************************************************************
clear;

%% Parameter settings
max_age = 1;                    % maximum age
age_interval = [0, max_age];    % age interval
space_interval = [-1, 1];       % spatial interval
discontinuity_age = 0.5;        % discontinuity in the age direction

pi0 = @(a) max_age - a;         % the survival rate

tau_array = 1;                  % boundary condition; tau_array = [0, 0.3, 1]
node_nummbers = 36;             % node array; node_nummbers = [6:2:60, 100]
d_array = 1e+0;                 % diffusion coefficient; d_array = 10.^(-10:8)

p = true;                       % plotting flag; set to false when looping over multiple data sets
if (length(tau_array) > 1 || length(node_nummbers) > 1 ...
                          || length(d_array) > 1)
    p = false;
end

%% β, μ1
% beta is independent of the spatial variable x, 
% mu1 is independent of the spatial variable x
% sigma = 0; rho = 0;

% beta is independent of the spatial variable x, 
% mu1 depends on the spatial variable x
% sigma = 0; rho = 1;

% beta depends on the spatial variable x, 
% mu1 depends on the spatial variable x
sigma = 1; rho = 1;

beta  = @(a, x) 8*(a >= discontinuity_age).*(2 - sigma*x.^2);
mu1   = @(a, x) 1./(1 + rho*x.^2);

% Folder for saving data
foldername = []; % ['.\example1_sigma', num2str(sigma), '_rho', num2str(rho), '\'];

%% Eigenvalues and eigenfunctions of the age-structured local diffusion model
for i = 1:length(tau_array)
    tau = tau_array(i);
    for j = 1:length(node_nummbers)
        N = node_nummbers(j);
        M = node_nummbers(j);
        for k = 1:length(d_array)
            d = d_array(k);
            [eigenvalues, eigenfunctions] = age_diffusion_model(...
                age_interval, space_interval, d, beta, pi0, mu1, ...
                tau, N, M, p, discontinuity_age, foldername);

            fprintf('--------------- d = %.12f completed! ---------------\n', d);
        end
    end
end

