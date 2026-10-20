%**************************************************************************
    %> Filename: example2.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-02-02 17:33:54
    %> Description: Numerical example of the eigenvalue problem for the
    %  age-structured local diffusion model with periodic boundaries.
%**************************************************************************
clear;

%% Parameter settings
max_age = 1;                    % maximum age
age_interval = [0, max_age];    % age interval
s = 0;                          % starting point in the spatial direction
discontinuity_age = 0.5;        % discontinuity in the age direction

pi0 = @(a) max_age - a;         % the survival rate

node_nummbers = 36;             % node array; node_nummbers = [6:2:60, 100]
r_array = 1;                    % values of the parameter r; r_array = 0.05:0.025:1
T_array = 1;                    % values of the period T; T_array = 10.^(-6:0.25:6);
d_array = 1e+0;                 % diffusion coefficient; d_array = 10.^(-10:8)

p = true;                       % plotting flag; set to false when looping over multiple data sets
if (length(node_nummbers) > 1 || length(r_array) > 1 ...
       || length(T_array) > 1 || length(d_array) > 1)
    p = false;
end

%% β, μ1
% beta is independent of the spatial variable x, 
% mu1 is independent of the spatial variable x
sigma = 0; rho = 0;

% beta is independent of the spatial variable x, 
% mu1 depends on the spatial variable x
% sigma = 0; rho = 1;

% beta depends on the spatial variable x, 
% mu1 depends on the spatial variable x
% sigma = 1; rho = 1;

beta  = @(a, x) 8*(a >= discontinuity_age).*(2 + sigma*sin(2*pi*x));
mu1   = @(a, x) exp(rho*cos(2*pi*x));

% Folder for saving data
foldername = []; % ['.\example2_sigma', num2str(sigma), '_rho', num2str(rho), '\'];

%% Eigenvalues and eigenfunctions of the age-structured local diffusion model
for i = 1:length(node_nummbers)
    N = node_nummbers(i);
    M = node_nummbers(i);
    for j = 1:length(r_array)
        r = r_array(j);
        for k = 1:length(T_array)
            T = T_array(k);
            for q = 1:length(d_array)
                d = d_array(q);
                [eigenvalues, ~] = age_diffusion_model_period(...
                    age_interval, s, r, T, d, beta, pi0, mu1, ...
                    N, M, p, discontinuity_age, foldername);
            end
        end

        fprintf('--------------- r = %.12f completed! ---------------\n', r);
    end
end

