clear all;clc;
m1=10;m2=350;
kw=500000;ks=10000;
b=500;

A=[0 1 0 0;
   -(kw+ks)/m1 -b/m1 ks/m1 b/m1;
   0 0 0 1;
   ks/m2 b/m2 -ks/m2 -b/m2];
B=[0 kw/m1 0 0]';
C=[1 0 0 0;
   0 0 1 0];
D=[0;0];

sys=ss(A,B,C,D);
H=tf(sys)       % transfer function of the suspension system

num=[7.143e04 1.429e06];
den=[1 51.43 5.103e04 7.143e04 1.429e06];
syms s;
syms omega real;                       % symbolic calculation 
f1 = poly2sym(num,s)/poly2sym(den,s);
f2 = subs(f1,s,1i*omega);
f2_real = simplify(real(f2));
f2_imag = simplify(imag(f2));          % substitution omega for s in the transfer function

H2w=f2*conj(f2);                     % magnitude of transfer function squared (complex numbers)

CR=0.336e-6;w=2;
om=(0:0.1:100);
Sxx=CR*(om).^-w;           % load 


x0 = 0:1000;              
FUN = matlabFunction(H2w);  
y = feval(FUN, x0);          % evalution of transfer function numerically

Syy=y.*Sxx;    % final formula
% 
% figure
% plot(om,y)
% title('Response of damping system');
% xlabel('omega (1/s)');
% ylabel('(H(omega))^2(omega)');
% xlim([0 3])
% 
% figure
% plot(om,Sxx)
% title('Power spectral density of loading');
% xlabel('omega (1/s)');
% ylabel('Sxx(omega)(m^2/Hz)');
% xlim([0 3])
% 
% figure
% plot(om,sqrt(Sxx)*1000)
% title('Amplitude of loading');
% xlabel('omega (1/s)');
% ylabel('Elevation amplitude scaling factor (mm/sqrt(Hz)) ');
% xlim([0 5])
% 
% figure
% plot(om,Syy)
% title('Power spectral density of response');
% xlabel('omega (1/s)');
% ylabel('Yxx(omega) ((m/s^2)^2/sqrt(Hz)) ');
% xlim([0 3])
% 
% figure
% plot(om,sqrt(Syy)*1000)
% title('Amplitude of response');
% xlabel('omega (1/s)');
% ylabel('Vertical vibration scaling factor ((mm/s^2)/sqrt(Hz)) ');
% xlim([0 3])

t = 0:0.01:10;

%scaling factors for random load and response
loadScale = sqrt(Sxx).*1000; %Load
RScale = sqrt(Syy).*1000; %Response

sumL = 0; sumR = 0;

for c1 = 1:50
    for c2 = 1:20

    
r = randi([2,200]); %random point in scaling
f = (r/10); %frequency at the random point
phase = rand*pi; %random phase of wave

    %Loading
aL = loadScale(r)*sqrt(f); %Amplitude at the frequency
wL = aL*sin(f*2*pi*t+phase); %Sine wave
sumL = sumL + wL; %Summing sine waves

    %Response
aR = RScale(r)*sqrt(f); %Amplitude at the frequency
wR = aR*sin(f*2*pi*t+phase); %Sine wave
sumR = sumR + wR; %Summing sine waves
    end
end

figure
subplot(2,1,1)
plot(t,sumL)
xlabel('time [s]');
ylabel('Road elevation [mm]');

subplot(2,1,2)
plot(t,sumR)
xlabel('time [s]');
ylabel('Car vertical displacement [mm]');

%for t = 1:3
    %m_R = mean(sumR)
    %%m_L = mean(sumL)
    %%std_R = std(sumR)
    %std_L = std(sumL)

%autocorr_R = xcorr(sumR);
%autocorr_L = xcorr(sumL);
%end

figure
subplot(2,1,1)
histogram(sumL)
xlabel('Road elevation [mm]');
ylabel('Probability');

subplot(2,1,2)
histogram(sumR)
xlabel('Car vertical displacement [mm]');
ylabel('Probability');

% polyfit(


%[p,S,mu] = polyfit(x,y,n) also returns mu, which is a two-element vector with centering and scaling values. mu(1) is mean(x), and mu(2) is std(x). Using these values, polyfit centers x at zero and scales it to have unit standard deviation,


