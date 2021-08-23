clear all;

%t = [0:60]; % Time interval
Fs = 100/10;
T = 1/Fs;
L = 1000;
t=(0:L-1)*T;

F = 1;
m = 1; c = 0.5; k = 1;
w_e = 1; E = 0.25;
w_d = 0.9682458;
a_1 = 1; a_2 = 0;
a = 1; w = 0.9; phi = 0;

%Curves
steady_state = a*sin((w*t)- phi);
transient = exp(1).^(-E*w_e*t).*(a_2*sin((w_d*t)));
u_exp = exp(1).^(-E*w_e*t).*(a_1*cos(w_d*t));

result = u_exp + steady_state;

% fft of stationary part
Fs = 100/10;
T = 1/Fs;
L = 100;
t1=(0:L-1)*T;
steady_state_s = a*sin((w*t1)- phi);
Y = fft(steady_state_s);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:L/2)/L;
plot(f,P1);
title('Frequency of stationary part with FFt');
xlabel('f (Hz)');
ylabel('|P1(f)|');
xlim([0,10]);ylim([0,1.1]);

% fft of entire signal
Fs2 = 100/10;
T2 = 1/Fs2;
L2 = 100;
t2=(0:L2-1)*T2;
signal = (exp(1).^(-E*w_e*t2)).*((a_1*cos(w_d*t2))+(a_2*sin((w_d*t2))))+a*sin((w*t2)- phi);
Y2 = fft(signal);
P22 = abs(Y2/L2);
P12 = P22(1:L2/2+1);
P12(2:end-1) = 2*P12(2:end-1);

f2 = Fs2*(0:L2/2)/L2;
plot(f2,P12);
title('Frequency of entire signal with FFt');
xlabel('f (Hz)');
ylabel('|P1(f)|');
%xlim([0,10]);ylim([0,1.2]);

figure
plot(t,steady_state)
hold on;
plot(t,u_exp)
hold on;
plot(t,result)
hold off
grid on;
xlabel('Time [s]');
ylabel('Amplitude');
title('Curves')
legend('Steady State','Transient','Stationary Part')
