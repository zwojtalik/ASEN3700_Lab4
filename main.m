clc
close all
clear

%% Printing figures
printfigs = 1; % Set to 1 for saving/printing figures

%% Getting variables
load("lab04_analysis_signal1.mat");

%% Prelab
fs = 2048;
t_pl = 0:1/fs:1-1/fs;
x_pl = cos(2*pi*t_pl*300);

N = length(x_pl);
xdft = fft(x_pl);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x_pl):fs/2;

plot(freq,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")
%% Problem 4.1 
% a
figure()
plot(t(1:100), x(1:100))
xlabel('Time (s)')
ylabel('X values')
title('4.1.a: X Vals. vs. Time (s)')
grid on

if printfigs == 1
    print('-dpng', '-r300')
end

% b
figure()
fs_b = Fs; % Arbitrarily picked, can change this
t_b = t(1:100); % Time vector
x_b = x(1:100); % Input signal
N_b = length(x_b(1:100)); % Done for convinience
dt = t(2) - t(1); % dt - used later

xdft_b = fft(x_b); % Finding the FFT of our singal
xdft_b = xdft_b(1:N_b/2+1); % Halving the signal to only keep positive components

psdx_b = (dt^2/t_b(1,end)) * (abs(xdft_b).^2); % Calculating the power spectrum density at input x
%psdx_b = (1/(fs_b*N)) * (abs(xdft_b).^2); % Second PSDX to confirm first one

psdx_b(2:end-1) = 2*psdx_b(2:end-1); % We halved the signal earlier and are now doubling the signal strength to account for it

freq_b = 0:fs_b/N_b:fs_b/2; % Creating the frequency vector for the plot

dB_rms_b = 10*log10(psdx_b); % Converting from Prms into dB

plot(freq_b,dB_rms_b)
grid on
title('4.1.b: Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency dB(Vrms^2/Hz)')

if printfigs == 1
    print('-dpng', '-r300')
end

% c
figure()
fs_c = Fs; % Arbitrarily picked, can change this
t_c = t; % Time vector
x_c = x; % Input signal
N_c = length(x_c); % Done for convinience
dt = t(2) - t(1); % dt - used later

xdft_c = fft(x_c); % Finding the FFT of our singal
xdft_c = xdft_c(1:N_c/2+1); % Halving the signal to only keep positive components

psdx_c = (dt^2/t_c(1,end)) * (abs(xdft_c).^2); % Calculating the power spectrum density at input x

psdx_c(2:end-1) = 2*psdx_c(2:end-1); % We halved the signal earlier and are now doubling the signal strength to account for it

freq_c = 0:fs_c/N_c:fs_c/2; % Creating the frequency vector for the plot

dB_rms_c = 10*log10(psdx_c); % Converting from Prms into dB

plot(freq_c,dB_rms_c)
grid on
title('4.1.c: Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency dB(Vrms^2/Hz)')

if printfigs == 1
    print('-dpng', '-r300')
end

%% 4.1 d
load("lab04_analysis_signal2.mat");
% d1
figure()
plot(t(1:100), x(1:100))
xlabel('Time (s)')
ylabel('X values')
title('4.1.d1: X Vals. vs. Time (s)')
grid on

if printfigs == 1
    print('-dpng', '-r300')
end

% d2
[freq,psdx] = fft_func(Fs, t(1:100), x(1:100));

figure()
plot(freq, 10*log10(psdx))
grid on
title('4.1.d2: Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency dB(Vrms^2/Hz)')

if printfigs == 1
    print('-dpng', '-r300')
end

% d3
[freq,psdx] = fft_func(Fs, t, x);

figure()
plot(freq, 10*log10(psdx))
grid on
title('4.1.d3: Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency dB(Vrms^2/Hz)')

if printfigs == 1
    print('-dpng', '-r300')
end

%% 4.1 e
load("lab04_analysis_signal3.mat");
% e1
figure()
plot(t(1:100), x(1:100))
xlabel('Time (s)')
ylabel('X values')
title('4.1.e1: X Vals. vs. Time (s)')
grid on

if printfigs == 1
    print('-dpng', '-r300')
end

% e2
[freq,psdx] = fft_func(Fs, t(1:100), x(1:100));

figure()
plot(freq, 10*log10(psdx))
grid on
title('4.1.e2: Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency dB(Vrms^2/Hz)')

if printfigs == 1
    print('-dpng', '-r300')
end

% e3
[freq,psdx] = fft_func(Fs, t, x);

figure()
plot(freq, 10*log10(psdx))
grid on
title('4.1.e3: Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency dB(Vrms^2/Hz)')

if printfigs == 1
    print('-dpng', '-r300')
end