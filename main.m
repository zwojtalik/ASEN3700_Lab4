clc
close all
clear

%% Getting variables
load("lab04_analysis_signal1.mat");

%% Problem 4.1 
% a
figure()
plot(t(1:100), x(1:100))
xlabel('Time (s)')
ylabel('X values')
title('X Vals. vs. Time (s)')
grid on