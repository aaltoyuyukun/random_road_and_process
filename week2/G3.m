%Time
clear all;

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 3000;             % Length of signal
t = (0:L-1)*T;        % Time vector


%Sinusoidal signals
%Amplitudes
A1=2;A2=1;A3=0.7;
%Frequencies
f1=4;f2=30;f3=100;
%Phase angles
phi1=pi/4;phi2=pi;phi3=0;

%Sin waves
SS1 = A1*sin(2*pi*f1*t+phi1);
SS2 = A2*sin(2*pi*f2*t+phi2);
SS3 = A3*sin(2*pi*f3*t+phi3);

%Sum of sin waves
Sum = SS1+SS2+SS3;

%Adding noise
Noise = 2*randn(size(t));
SumN = Sum + Noise;

figure
subplot(2,1,1);
hold on
plot(t,SS1)
plot(t,SS2)
plot(t,SS3)
hold off
grid on;
ylabel('amplitude');
xlabel('Time [s]');
title('Every sinusoidal waves')
legend('SS1','SS2','SS3')

subplot(2,1,2)
plot(t,Noise)
grid on;
ylabel('amplitude');
xlabel('Time [s]');
title('Signal corruction with random noise')


figure
plot(t,SumN)
grid on;

ylabel('amplitude');
xlabel('Time [s]');
title('Sum of sinusoidal waves')

%Computing fast fourier transform

%FFT with noise
Y = fft(SumN);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure
plot(f,P1)
grid on;
title('FFT of wave signal with noise')
xlabel('f (Hz)')
ylabel('Amplitude')

%Narrow band process

%Sinusoidal signals
%Amplitudes
An1=3;An2=2;An3=1;
%Frequencies
fn1=5;fn2=5.5;fn3=6;
%Phase angles
phin1=pi/3;phin2=pi/2;phin3=0;

%Sin waves
SSn1 = An1*sin(2*pi*fn1*t+phin1);
SSn2 = An2*sin(2*pi*fn2*t+phin2);
SSn3 = An3*sin(2*pi*fn3*t+phin3);

%Sum of sin waves
narrowSum = SSn1+SSn2+SSn3;

figure
subplot(2,1,1);
hold on
plot(t,SSn1)
plot(t,SSn2)
plot(t,SSn3)
hold off
grid on;
ylabel('amplitude');
xlabel('Time [s]');
title('Every sinusoidal waves in narrow band')
legend('SSn1','SSn2','SSn3')

subplot(2,1,2)
plot(t,narrowSum)
grid on;
ylabel('amplitude');
xlabel('Time [s]');
title('Sum of sinusoidal signals in Narrow band ')

%FFT on Narrow band process
Yn = fft(narrowSum);
Pn2 = abs(Yn/L);
Pn1 = Pn2(1:L/2+1);
Pn1(2:end-1) = 2*Pn1(2:end-1);

figure
plot(f,Pn1)
grid on;
title('FFT of wave signal with narrow band signal')
xlabel('f (Hz)')
ylabel('Amplitude')

%Broad band process
broadSum = zeros(size(t));
%Sin waves
for c = 1:5
    broadSum = broadSum + 10*rand*sin(2*pi*rand*100*t+rand*2*pi);
end

%Sum of sin waves
figure
plot(t,broadSum)
grid on;
ylabel('amplitude');
xlabel('Time [s]');
title('Sum of sinusoidal signals in broad band ')

%FFT on Narrow band process
Yb = fft(broadSum);
Pb2 = abs(Yb/L);
Pb1 = Pb2(1:L/2+1);
Pb1(2:end-1) = 2*Pb1(2:end-1);

figure
plot(f,Pb1)
grid on;
title('FFT of wave signal with broad band signal')
xlabel('f (Hz)')
ylabel('Amplitude')

%Using hand window to improve the signal

w = hann(L);

wnarrowSum = zeros(size(w));
wbroadSum = zeros(size(w));

for c1 = 1:L
    WnarrowSum(c1) = narrowSum(c1)*w(c1);
    WbroadSum(c1) = broadSum(c1)*w(c1);
end

figure
plot(w)
grid on;
title('Hann windowing')
xlabel('Samples')
ylabel('Amplitude')

figure
plot(t,WnarrowSum)
grid on;
title('Sum of narrow band waves with windowing')
xlabel('Time [t]')
ylabel('Amplitude')

figure
plot(t,WbroadSum)
grid on;
title('Sum of broad band waves with windowing')
xlabel('Time [t]')
ylabel('Amplitude')

wYn = fft(WnarrowSum);
wPn2 = abs(wYn/L);
wPn1 = wPn2(1:L/2+1);
wPn1(2:end-1) = 2*wPn1(2:end-1);

wYb = fft(WbroadSum);
wPb2 = abs(wYb/L);
wPb1 = wPb2(1:L/2+1);
wPb1(2:end-1) = 2*wPb1(2:end-1);

figure
plot(f,wPn1)
grid on;
title('FFT on Sum of narrow band waves with windowing')
xlabel('f [Hz]')
ylabel('Amplitude')

figure
plot(f,wPb1)
grid on;
title('FFT on Sum of broad band waves with windowing')
xlabel('f [Hz]')
ylabel('Amplitude')