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
% sampling frequency,linspace(x1,x2,N)
om=linspace(2*pi/90.9,2*pi/0.3,1001); 
Sxx=CR*(om).^-w; % load 

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
% % 
% figure
% plot(om,sqrt(Sxx)*1000)
% title('Amplitude of loading');
% xlabel('omega (1/m)');
% ylabel('Elevation amplitude scaling factor (mm/sqrt(Hz)) ');
% xlim([0 3])

% figure
% plot(om,Syy)
% title('Power spectral density of response');
% xlabel('omega (1/s)');
% ylabel('Yxx(omega) ((m/s^2)^2/sqrt(Hz)) ');
% xlim([0 3])
% % 
% figure
% plot(om,sqrt(Syy)*1000)
% title('Amplitude of response');
% xlabel('omega (1/m)');
% ylabel('Car vertical displacement scaling factor (mm/sqrt(Hz)) ');
% xlim([0 3])
%%%%%%%%%%%%%%%%

t = 0:0.01:50;

%scaling factors for random load and response
loadScale = sqrt(Sxx).*1000; %Load
RScale = sqrt(Syy).*1000; %Response

sumL = 0; sumR = 0;

samples = 1001; % for fitting
for c1 = 1:samples
    
r = randi([2,1001]); %random point in scaling
f = (r/354); %frequency at the random point
phase = rand*2*pi; %random phase of wave

    %Loading
aL = loadScale(r)*sqrt(f); %Amplitude at the frequency
wL = aL*sin(f*2*pi*t+phase); %Sine wave
sumL = sumL + wL; %Summing sine waves

    %Response
aR = RScale(r)*sqrt(f); %Amplitude at the frequency
wR = aR*sin(f*2*pi*t+phase); %Sine wave
sumR = sumR + wR; %Summing sine waves
end

% sumL = sumL/samples;
% sumR = sumR/samples;


figure
subplot(2,1,1)
plot(t,sumL)
xlabel('Travel [m]');
ylabel('Road elevation [mm]');

subplot(2,1,2)
plot(t,sumR)
xlabel('Travel [m]');
ylabel('Car vertical displacement [mm]');

figure
histogram(sumL,'Normalization','pdf')
% figure
% histogram(sumR,'Normalization','pdf')




