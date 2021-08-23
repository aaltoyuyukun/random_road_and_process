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
om=linspace(2*pi/90.9,2*pi/0.3,1001);
Sxx=CR*(om).^-w;           % load 

x0 = 0:1000;              
FUN = matlabFunction(H2w);  
y = feval(FUN, x0);          % evalution of transfer function numerically

Syy=y.*Sxx;    % final formula

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
%%%%%%%%%%%%%%%

t = 0:0.001:100;

%scaling factors for random load and response
loadScale = sqrt(Sxx).*1000; %Load
RScale = sqrt(Syy).*1000; %Response

sumL = 0; sumR = 0;

samples = 1001;
for c1 = 1:samples
    
r = randi([2,1001]); %random point in scaling
f = (r/354); %frequency at the random point
phase = rand*2*pi; %random phase of wave

    %Loading
aL = loadScale(r)*sqrt(f); %Amplitude at the frequency
wL = aL*sin(f*2*pi*t+phase); %Sine wave
sumL = sumL + wL; %Summing sine waves
end

%Response
Response = lsim(sys,sumL,t); 
% sumL is artificial road wave
% Response is the response from the system after state function

% figure
% subplot(2,1,1)
% plot(t,sumL)
% xlabel('Travel [m]');
% ylabel('Road elevation [mm]');
% 
% subplot(2,1,2)
% plot(t,Response(:,2))
% xlabel('Travel [m]');
% ylabel('Car vertical displacement [mm]');

% figure
% histogram(sumL,'Normalization','pdf')
% title('Continuous distribution - Loading')
% 
% figure
% histogram(Response(:,2),'Normalization','pdf')
% title('Continuous distribution - Response')

%%%%% FITTING THE CURVE TO HISTOGRAM

max = max(sumL); min = min(sumL);            %finding min and max value 
gap = 0.1;                                   %define gap between levels
num_lev = (max-min)/gap;                     %define the number of levels

% calculating probability exceed certain road elevation, after iteration level is lowered by the gap
for i = 1:num_lev
layer(i) = min + (i * gap); % 逐层递进
count(i) = sum(sumL > layer(i)); % 挑出计数sumL大于第i层值的数量
cumula_prob(i) = count(i)/numel(t); % 计算概率
prob_vector(i) = min + (i * gap); 
end

num3 = numel(prob_vector) - 1;     % density function points based on random data
k = zeros(1,num3);
for a = 1:num3
k(a) = -(cumula_prob(a+1)-cumula_prob(a))/(prob_vector(a+1)-prob_vector(a));
end
k(num3+1) = 0;

% figure
% histogram(sumL,'Normalization','pdf')
% title('Continuous distribution - Loading')
% hold on;

x = 1:numel(k);                   % fit curve on calculated distribution points
x2 = polyfit(prob_vector,k,9);
y3 = polyval(x2,prob_vector);
standart_deviation_load = std(sumL);
% pdf('Normal',随机量,mean_value=0,标准差))
y2 = pdf('Normal',prob_vector,0,standart_deviation_load);
plot(prob_vector,y2,'r','LineWidth',1.2);
hold on;
plot(prob_vector,y3,'k','LineWidth',1.2);
legend('Histogram','Normal','Polyfit');

% figure
% histogram(Response(:,2),'Normalization','pdf');
% title('Continuous distribution - Response');
% hold on;

standart_deviation_response = std(Response(:,2));
y2 = pdf('Normal',prob_vector,0,standart_deviation_response);
plot(prob_vector,y2,'r','LineWidth',1.2);
legend('Histogram','Normal');

% [c,hist,edges,rmm,idx] = rainflow(Response(:,2),1/t);
% T = array2table(c,'VariableNames',{'Count','Range','Mean','Start','End'});
% 
% figure;
% histogram('BinEdges','edges','BinCounts',sum(hist,2));
% title('Cycle counts as a function of stress');
% xlabel('Stress Range');
% ylabel('Cycle Counts');
% 
% figure;
% rainflow(Response,1/t);





