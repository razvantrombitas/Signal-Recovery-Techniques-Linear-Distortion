% Clear workspace, command window, and close all figures
close all
clear all
clc

%A. Read the audio file

[s0, fs] = audioread('signal.wav');

% plot(s0,'b','LineWidth',1.5);
% title('Original Signal');
% ylabel('s(n)');
% xlabel('(n)');
% figure

%B. Normalize the signal

s1 = s0/max(abs(s0));

% plot(s1,'b','LineWidth',1.5);
% title('Normalized Signal');
% ylabel('s1(n)');
% xlabel('n');


%C.
fe = 16000; 
s2 = resample(s1,160,441); 

% plot(s2,'r','LineWidth',1);
% title('Resampled Signal');
% ylabel('s2(n)');
% xlabel('(n)');

% s_interpolation = interp(s1, 160);
% s_decimation = decimate(s_interpolation, 441);

%D. Apply a Butterworth filter
fc = 7000; 
% Pass through a low-pass filter (a Butterworth filter is chosen because its frequency response is smooth - without ripples - both in the passband and the stopband)

[b,a] = butter(10,2*fc/fe);
[H,w] = freqz(b,a,[],fe);

s = filter(b,a,s2);

% s = s/max(abs(s)); % normalize for the use of the audiowrite function without warning
audiowrite('filter_d.wav',s,fe);

%E. MATLAB implementation of the digital distortion filter using B and A
%coefficients
B = [0.4915, 0.1826, 1.4321, 0.9146, 2.1694, 1.9138, 2.2428, 2.4305, 1.7964, 2.0088, 1.1784, 1.0973, 0.6063, 0.3872, 0.2106, 0.0740, 0.0379];
A = [1.0000, -0.6452, 2.1976, -0.3915, 1.9283, 1.0569, 0.4561, 2.0150, -0.1674, 1.3202, 0.1298, 0.3213, 0.2968, -0.0493, 0.1237, -0.0145, 0.0093];
x = filter(B,A,s);

% x = x/max(abs(x)); % normalize for the use of the audiowrite function without warning
audiowrite('filter_e.wav',x,fe);

%F. Plot the frequency response of the original signal and the filtered signal
%frequency
frecv_s = (1:length(s))*fe/length(s);
S = fft(s);
plot(frecv_s,abs(S),'b','LineWidth',1);
title('Signal s');
xlabel('Frequency (Hz)')
ylabel('Magnitude response of signal s');

hold on

frecv_x = (1:length(x))*fe/length(x);
X = fft(x);
plot(frecv_x,abs(X),'b','LineWidth',1);
title('Distorted signal x');
xlabel('Frequency (Hz)')
ylabel('Magnitude response of signal x');

figure

% Plot the signal s
plot(s,'b','LineWidth',1);
title('Signal s');
ylabel('s(n)');
xlabel('n');

figure;

% Plot the distorted signal x
plot(x,'r','LineWidth',1);
title('Distorted signal x');
ylabel('x(n)');
xlabel('n');

figure;

% Plot the squared error between s and x
plot((s-x).^2)
title('Squared error between signals s and x');
ylabel('SE');
xlabel('n');

figure;

% Generate the frequency response of the distortion filter
[H_filtru,w_filtru] = freqz(B,A,[],fe);
plot(w_filtru,abs(H_filtru))
title('Distortion filter frequency response');
ylabel('Magnitude of frequency response of distortion filter');
xlabel('Frequency (Hz)')

figure;

% Generate white noise
t_n = (1:10000);
s_aux = randn(size(t_n));
x_aux = filter(B,A,s_aux);

% s_aux -> desired signal passed through the distortion filter
% x_aux -> output signal from the distortion filter
% y_aux -> "clean" signal at the output of the adaptive filter (similar to the input signal s_aux)

% Set parameters for the NLMS adaptive filter
mu = 0.2; %mu - step size (0 < mu < 1; controls stability and convergence speed of the algorithm)
sigma = 1; %estimated power of the input signal
alpha = 0; %correction factor (0 < alpha << 1);
px = 0; %vector used by the function to store previous values of the input signal
L = 15; %L - order of the adaptive filter (FIR)
b = zeros(1, L+1); %b - vector of filter coefficients (L + 1 coefficients)

% Apply the NLMS adaptive filter to the distorted signal x_aux and the desired signal s_aux
[y_aux,b,px] = nlms(x_aux,s_aux,b,mu,sigma,alpha,px);

%H.
%frequency
frecv_s_aux = (1:length(s_aux))*fe/length(s_aux);
S_aux = fft(s_aux);
plot(frecv_s_aux,abs(S_aux))
title('Frequency response of signal s_a_u_x');
ylabel('Magnitude of signal s_a_u_x response');
xlabel('Frequency (Hz)')

figure

frecv_x_aux = (1:length(x_aux))*fe/length(x_aux);
X_aux = fft(x_aux);
plot(frecv_x_aux,abs(X_aux))
title('Frequency response of signal x_a_u_x');
ylabel('Magnitude of signal x_a_u_x response');
xlabel('Frequency (Hz)')

figure

frecv_y_aux = (1:length(y_aux))*fe/length(y_aux);
Y_aux = fft(y_aux);
plot(frecv_y_aux,abs(Y_aux))
title('Frequency response of signal y_a_u_x');
ylabel('Magnitude of signal y_a_u_x response');
xlabel('Frequency (Hz)')

figure

[H_filtru_adaptiv,w_filtru_adaptiv] = freqz(b,1,[],fe);
plot(w_filtru_adaptiv,abs(H_filtru_adaptiv))
title('Adaptive filter frequency response');
ylabel('Magnitude of adaptive filter response in frequency');
xlabel('Frequency (Hz)')

figure

plot(s_aux)
title('Time domain plot of signal s_a_u_x');
ylabel('s_a_u_x(n)');
xlabel('(n)')

figure

plot(x_aux)
title('Time domain plot of signal x_a_u_x');
ylabel('x_a_u_x(n)');
xlabel('(n)')

figure

plot(y_aux)
title('Time domain plot of signal y_a_u_x');
ylabel('y_a_u_x(n)');
xlabel('(n)')

figure

%Obs. The two frequency responses multiplied point by point should give 1.

plot(w_filtru,ones(size(w_filtru)),'r','LineWidth',1)
hold on

produs = abs(H_filtru).*abs(H_filtru_adaptiv);
plot(w_filtru,produs)

legend ('Ideal','Real');
title('Product of frequency responses of adaptive filter and distortion filter');
ylabel('Product');
xlabel('Frequency (Hz)')

figure

% I.
y = filter(b, 1, x);
audiowrite('filter_i.wav',y,fe);

% J.
% frequency
frecv_s = (1:length(s))*fe/length(s);
S = fft(s);
plot(frecv_s,abs(S),'b','LineWidth',1);
title('Frequency response of signal s');
ylabel('Magnitude of response of signal s');
xlabel('Frequency (Hz)')

hold on

frecv_y = (1:length(y))*fe/length(y);
Y = fft(y);
plot(frecv_y,abs(Y),'r','LineWidth',1)
title('Frequency response of signal y');
ylabel('Magnitude of response of signal y');
xlabel('Frequency (Hz)')

legend('signal s','signal y');

figure

frecv_x = (1:length(x))*fe/length(x);
X = fft(x);
plot(frecv_x,abs(X))
title('Frequency response of signal x');
ylabel('Magnitude of response of signal x');
xlabel('Frequency (Hz)')

figure

% time
plot(s)
title('Signal s in time');
ylabel('s(n)');
xlabel('(n)')

hold on

plot(y)
title('Signal y in time');
ylabel('y(n)');
xlabel(' (n)')

legend('signal s','signal y');

figure

plot(x)
title('Signal x in time');
ylabel('x(n)');
xlabel('(n)')

% Observation: The degree of similarity between signals is high after optimizing the adaptive filter

% Initial values:
% mu = 0.08;
% sigma = 1;
% alpha = 0;
% px = 0;
% L = 20;
% b = zeros(1, L+1);

% Values after optimization:
% mu = 0.2;
% sigma = 1;
% alpha = 0;
% px = 0;
% L = 15;
% b = zeros(1, L+1);

