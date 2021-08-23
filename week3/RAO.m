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
H=tf(sys);

num=[7.143e04 1.429e06];
den=[1 51.43 5.103e04 7.143e04 1.429e06];
syms s;
syms omega real;
f1 = poly2sym(num,s)/poly2sym(den,s); % 返回多项式系数
f2 = subs(f1,s,1i*omega); % subs(a+b,a,4) = b + 4，替换频率域到复频域
f2_real = simplify(real(f2));
f2_imag = simplify(imag(f2));

H2w=f2*conj(f2); % 计算共轭复数 Conjugate complex number

CR=0.336e-6;w=2;
om=(0:0.1:100);
Sxx=CR*(om).^-w;

x0 = 0:1000;
FUN = matlabFunction(H2w);  % 指定用于计算的函数
y = feval(FUN,x0);
% 得到1-1000个值，将这些值往FUN函数里面代入，y即为频率域图像

Syy=y.*Sxx;
plot(om,y);
title('Response of damping system');
xlabel('omega (1/s)');
ylabel('(H(omega))^2(omega)');
xlim([0 1])