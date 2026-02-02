function [freq,psdx] = fft_func(Fs, t, x)
%FFT Takes inputs of:
%   Fs = frequency
%   t = Time vector
%   x = Input signal
% 
% And returns outputs of:
%   freq = Frequency array for plot
%   psdx = Power Spectrum Density 

N = length(x); % Done for convinience
dt = t(2) - t(1); % dt - used later

xdft = fft(x); % Finding the FFT of our singal
xdft = xdft(1:N/2+1); % Halving the signal to only keep positive components

psdx = (dt^2/t(1,end)) * (abs(xdft).^2); % Calculating the power spectrum density at input x

psdx(2:end-1) = 2*psdx(2:end-1); % We halved the signal earlier and are now doubling the signal strength to account for it

freq = 0:Fs/N:Fs/2; % Creating the frequency vector for the plot

end