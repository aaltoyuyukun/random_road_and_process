% part3 freak event
Fs = 1000;
T = 1/Fs;
L = 1000;
t = (0:L-1) * 10 * T;

A = 2 * ones(1,10);
w = rand * [5.1:0.2*2:7];
phi_freak = [(-pi+pi/10):rand * pi/2.5:pi];
S_freak = 0;
sum_S_freak = 0;
sum_X_freak = 0;
figure;
subplot(2,1,1);
for j=1:length(w)
    S_freak= A(j) * sin(2 * pi * w(j) * t + rand^2 * phi_freak(j)); %sum of signal without disturbance
    % X = S + 2 * rand(size(t)); %sum of signal with disturbance
    plot(t,S_freak);
    grid on;
    if j == length(-1)
        title('multiple sinusoidal signals including freak event in time domain');
        xlabel('s');
        ylabel('amplitude');
        hold off;
    end
    hold on;
    sum_S_freak = sum_S_freak + S_freak;
    % sum_X = sum_X + X;
end
subplot(2,1,2);
plot(t,sum_S_freak);
title('freak event');
xlabel('s');
ylabel('amplitude');
grid on;
