function [ X, freq ] = fft_func(x, Nf)
% ME 621                                                        Fall 2017
% compute FFT of signal and frequency domain (Hz)

% time & frequency domains
N    = size(x,2);                   % # data points
dt   = 1/(2*Nf);                    % time step (s)
T    = N*dt;                        % time period (s)
df   = 1/T;                         % frequency spacing (Hz)
freq = 0:df:Nf;                     % frequency domain (Hz)
X    = fft(x,[],2);                 % double-sided Fourier transform
X    = (2/N)*X(:,1:length(freq),:); % scaled single-sided