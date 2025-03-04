% 置信度及非线性函数应用:

clc;clear;close all;
addpath('Fncs\')

% 接收信号 y

y=rand(1000,1);


gamma=reliability(y);

a=5;
b=0.2;
f=reliability_nonlinear(gamma,a,b);