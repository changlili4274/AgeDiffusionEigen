%**************************************************************************
    %> Filename: plot_error_example2.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-02-06 17:27:15
    %> Description: Visualize the variation of the principal eigenvalue
    %  errors of Example 2 with the number of grid points.
%**************************************************************************
clear;

%% Parameter settings
temp_array = 6:2:60;
node_nummbers = [temp_array, 100];
node_nummbers2 = [60*ones(size(temp_array)), 100];
r = 1;
T = 1;
d = 1;

% Output figure folder
outfolder = './figures_example2/';
if ~isfolder(outfolder)
    mkdir(outfolder);
end

%% Load data
sigma = 1; rho = 1;
infolder = ['.\example2_sigma', num2str(sigma), '_rho', num2str(rho), '\'];

lambda1_array = zeros(length(node_nummbers), 1);
lambda3_array = zeros(length(node_nummbers), 1);
lambda5_array = zeros(length(node_nummbers), 1);

for j = 1:length(node_nummbers)
    N = node_nummbers2(j);
    M = node_nummbers(j);

    filename = [infolder, 'N', num2str(N), '_M', num2str(M), ...
        '_r', num2str(r), '_T', num2str(T), '_d', num2str(d), '.mat'];
    load(filename, 'lambda1', 'lambda3', 'lambda5');

    lambda1_array(j) = lambda1;
    lambda3_array(j) = lambda3;
    lambda5_array(j) = lambda5;
end

% Compute the errors
lambda1_error = abs(lambda1_array(1:end-1) - lambda1_array(end));
lambda3_error = abs(lambda3_array(1:end-1) - lambda3_array(end));
lambda5_error = abs(lambda5_array(1:end-1) - lambda5_array(end));

%% Visualization
x = 0:0.5:60;
rate = exp(- 2.3*x);

figure(1);

semilogy(node_nummbers(1:end-1), lambda1_error(1:end), ...
    'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 14);
hold on;
semilogy(node_nummbers(1:end-1), lambda3_error(1:end), ...
    'g-s', 'LineWidth', 2, 'MarkerFaceColor', 'g', 'MarkerSize', 12);
hold on;
semilogy(node_nummbers(1:end-1), lambda5_error(1:end), ...
    'b-d', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 10);
hold on;
semilogy(x, rate, 'k--', 'LineWidth', 2);

hold off;
grid on;

% xlabel('\it{N = M}');
xlabel('\it{M (N = 60)}');
% xlabel('\it{N (M = 60)}');

ylabel('\it{error}');
axis([0, 60, 1e-13, 1e+2])
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
set(get(gca,'Children'), 'linewidth', 2.5);
set(get(gca,'YLabel'), 'FontSize', 25, 'FontName', 'Times New Roman');
set(get(gca,'XLabel'), 'FontSize', 25, 'FontName', 'Times New Roman');
legend({'\lambda_1', '\lambda_3', '\lambda_5', '\it{e^{-2.3M}}'}, ...
        'Location', 'northeast');

% filename = [outfolder, 'convergence_sigma', num2str(sigma), ...
%     '_rho', num2str(rho)];

filename = [outfolder, 'convergence_sigma', num2str(sigma), ...
    '_rho', num2str(rho), '_N60'];

% filename = [outfolder, 'convergence_sigma', num2str(sigma), ...
%     '_rho', num2str(rho), '_M60'];
    
print(gcf, filename, '-djpeg', '-r600');
saveas(gcf, filename, 'fig');

