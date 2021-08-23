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
H=tf(sys);       % transfer function of the suspension system

num=[7.143e04 1.429e06];
den=[1 51.43 5.103e04 7.143e04 1.429e06];
syms s;
syms omega real;                       % symbolic calculation 
f1 = poly2sym(num,s)/poly2sym(den,s);
f2 = subs(f1,s,1i*omega);
% f2_real = simplify(real(f2));
% f2_imag = simplify(imag(f2));          % substitution omega for s in the transfer function

H2w=f2*conj(f2);                     % magnitude of transfer function squared (complex numbers)

CR=0.336e-6;w=2;
om=(0:0.1:100);
Sxx=CR*(om).^-w;           % load response

x0 = 0:1000;              
FUN = matlabFunction(H2w);  
y = feval(FUN, x0);          % evalution of transfer function numerically

Syy=y.*Sxx;    % final formula
plot(om,y)
title('Response of damping system');
xlabel('omega (1/s)');
ylabel('(H(omega))^2(omega)');
xlim([0 1])
