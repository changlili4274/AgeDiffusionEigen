%**************************************************************************
    %> Filename: fourier_difmat.m
    %> Author: Lili Chang
    %> Mail: changll@amss.ac.cn
    %> Created Time: 2026-03-18 16:14:19
    %> Description: Compute the first-, second-, and third-order
    %  differentiation matrices for the Fourier spectral method.
%**************************************************************************
function [x, D1, D2, D3] = difmat_fourier(a, b, N)
%   Input:
%       a, b - interval endpoints; periodic boundary conditions; grid points on [a, b)
%       N    - number of grid nodes (positive integer)
%   Output:
%       x    - column vector of length N containing the grid nodes
%       D1, D2, D3 - N-by-N differentiation matrices (first, second, and third order)
%
%   Note: For even N, the wavenumber corresponding to the Nyquist frequency
%         is set to zero in the first- and third-order derivatives to avoid
%         complex results; the second-order derivative retains that frequency
%         (because its square is real).

%%
if nargin < 3 || a >= b
    error('a must be less than b.');
end

if N < 2 || fix(N) ~= N
    error('N must be an integer >= 2 (number of nodes).');
end

L = b - a;              % interval length
h = L / N;              % grid spacing
x = a + (0:N-1)' * h;   % grid nodes (M-by-1)
    
% Unit impulse vector e0 = [1; 0; ...; 0]
e0 = zeros(N, 1);
e0(1) = 1;
    
%% Wavenumber definitions
% Wavenumbers for the first- and third-order derivatives (odd operators; 
% Nyquist frequency set to zero when M is even)
if mod(N, 2) == 0   % even
    % positive frequencies: 0, 1, ..., M/2-1
    % negative frequencies: -M/2+1, ..., -1
    % Nyquist frequency (M/2) set to zero
    k1 = 2*pi/L * [0:N/2-1, 0, -N/2+1:-1]';
else                % odd
    k1 = 2*pi/L * [0:(N-1)/2, -(N-1)/2:-1]';
end
    
% Wavenumbers for the second-order derivative (even operator; Nyquist 
% frequency retained when M is even)
if mod(N, 2) == 0   % even
    % positive frequencies: 0, 1, ..., M/2   (including M/2)
    % negative frequencies: -M/2+1, ..., -1
    k2 = 2*pi/L * [0:N/2, -N/2+1:-1]';
else                % odd
    k2 = k1;
end
    
%% Compute the first column
U0 = fft(e0);
    
% First-order derivative operator: D1 = ifft( i*k .* fft(u) )
c1 = real(ifft(1i * k1 .* U0));
    
% Second-order derivative operator: D2 = ifft( -k^2 .* fft(u) )
c2 = real(ifft(- (k2.^2) .* U0));
    
% Third-order derivative operator: D3 = ifft( -i*k^3 .* fft(u) )
c3 = real(ifft(-1i * (k1.^3) .* U0));
    
%% Construct circulant matrices
% The first row of a circulant matrix is obtained by a circular reversal 
% of the first column
first_row1 = [c1(1); c1(end:-1:2)]';
D1 = toeplitz(c1, first_row1);
    
first_row2 = [c2(1); c2(end:-1:2)]';
D2 = toeplitz(c2, first_row2);
    
first_row3 = [c3(1); c3(end:-1:2)]';
D3 = toeplitz(c3, first_row3);

end
