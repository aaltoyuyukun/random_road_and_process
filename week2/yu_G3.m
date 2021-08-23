%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total 3 parts
% part1 
% create multiple sinusodial functions and sum them
% perform fft and assess correctness
clear all;
Fs = 1000;
T = 1/Fs;
L = 1000;
t = (0:L-1) * T;

A = [1,2,3];
w = [30,60,120];
S = 0;
sum_S = 0;
sum_X = 0;
figure;
subplot(2,1,1);
for i=1:length(w)
    S= A(i) * sin(2 * pi * w(i) * t); %sum of signal without disturbance
    X = S + 2 * rand(size(t)); %sum of signal with disturbance
    plot(t,X);
    if i == length(-1)
        title('multiple sinusoidal signals in time domain');
        xlabel('s');
        ylabel('amplitude');
        hold off;
    end
    hold on;
    sum_S = sum_S + S;
    sum_X = sum_X + X;
end

subplot(2,1,2);
plot(t,sum_X,'r-');
title('sum of multiple sinusoidal signals in time domain');
xlabel('s');
ylabel('amplitude');

figure;
subplot(2,1,1);
Y = fft(sum_X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2 * P1(2:end-1); % 500 points

f = Fs * (0:L/2)/L;
plot(f,P1);
title('fft related to multiple sinusoidal signals with disturbance');
xlabel('f(Hz)');
ylabel('|P1(f)|');

subplot(2,1,2);
Y = fft(sum_S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2 * P1(2:end-1); % 500 points

f = Fs * (0:L/2)/L;
plot(f,P1);
title('fft related to multiple sinusoidal signals without disturbance');
xlabel('f(Hz)');
ylabel('|P1(f)|');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part2
% Show both narrow and broad banded processes in time and frequency domain
Fs = 1000;
T = 1/Fs;
L = 1000;
t = (0:L-1) * T;
tN = (0:L-1) * 3 * T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [1,2,3];
w = [3.5,4,4.5];
phi = [0,pi/3,pi/6,];
S = 0;
sum_S = 0;
sum_X = 0;
for i=1:length(w)
    % sum of signal without disturbance
    S = A(i) * sin(2 * pi * w(i) * tN + phi(i));
    % X = S + 2 * rand(size(t));
    sum_X = sum_X + S;
    % sum_S = sum_S + X;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nband = sum_X;
figure;
subplot(3,1,1);
plot(tN,Nband);
title('narrow banded process in time domain');
xlabel('s');
ylabel('amplitude');
grid on;
hold on;

YN = fft(Nband);
P2N = abs(YN/L);
P1N = P2N(1:L/2+1);
P1N(2:end-1) = 2 * P1N(2:end-1); % 500 points

subplot(3,1,2);
f = Fs * (0:L/2)/L;
plot(f,P1N);
hold off;
title('fft related to multiple sinusoidal signals with disturbance');
xlabel('f(Hz)');
ylabel('|P1(f)|');
grid on;

tNW = (0:L-1)*T;
WindowFN = hann(length(tNW));
NbandW = Nband .* WindowFN';

YNW = fft(NbandW);
P2NW = abs(YNW/L);
P1NW = P2NW(1:L/2+1);
P1NW(2:end-1) = 2 * P1NW(2:end-1); % 500 points

fNW = Fs * (0:L/2)/L;
subplot(3,1,3);
plot(fNW,P1NW);
title('fft related to multiple sinusoidal signals with WF');
xlabel('f(Hz)');
ylabel('|P1(f)|');
grid on;

figure;
plot(tNW,NbandW);
title('Narrow banded process with WindowFN');
xlabel('s');
ylabel('amplitude');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Broad band process @Vatii Lin
Fs = 1000;
T = 1/Fs;
L = 1000;
tB = (0:L-1) * T;
broadSum = zeros(size(tB));
% Sin waves
for c = 1:5
    broadSum = broadSum + 10*rand*sin(2*pi*rand*100*tB+rand*2*pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bband = broadSum;
figure;
subplot(3,1,1);
plot(tB,Bband);
title('borrow banded process in time domain');
xlabel('s');
ylabel('amplitude');
grid on;

YB = fft(Bband);
P2B = abs(YB/L);
P1B = P2B(1:L/2+1);
P1B(2:end-1) = 2 * P1B(2:end-1); % 500 points

f = Fs * (0:L/2)/L;
subplot(3,1,2);
plot(f,P1B);
title('fft related to multiple sinusoidal signals with broad banded process');
xlabel('f(Hz)');
ylabel('|P1(f)|');
grid on;

tBW = (0:L-1) * T;
WindowFB = hann(length(tBW));
BbandW = Bband .* WindowFB';

YNB = fft(BbandW);
P2NB = abs(YNB/L);
P1NB = P2NB(1:L/2+1);
P1NB(2:end-1) = 2 * P1NB(2:end-1); % 500 points

fWB = Fs * (0:L/2)/L;
subplot(3,1,3);
plot(fWB,P1NB);
title('fft related to multiple sinusoidal signals with WF');
xlabel('f(Hz)');
ylabel('|P1(f)|');
grid on;

figure;
plot(tBW,BbandW);
title('broad banded process with WindowFB');
xlabel('s');
ylabel('amplitude');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


