% 泰勒展开 数值结果对比
clc;clear;
addpath("Fncs\")

fs=20e9;
Ta=1/fs;

N=1e4;
% 
% 频率
f1=5e9;
f2=2e9;
f3=f1-f2;

% 载波强度和dither强度
a=0.5;
b=4;

% 信号和dither
z1=Creat_dither(fs,f1,N);
z2=Creat_dither(fs,f2,N*f2/f1);
z3=Creat_dither(fs,f3,N*f3/f1);


% 合并项
cos_term=2*b*z1+2*b*a*z2+2*a*z3;
const=b.^2+1+a.^2;
I=const+cos_term;

% ln项
x=log(sqrt(I));
y=1/2*(log(const)+cos_term./const);

% 误差
error=x-y;

% 检测是否满足展开条件
figure;
plot(cos_term./const)