% 泰勒展开 数值结果对比
clc;clear;

fs=10e9;
Ta=1/fs;

N=1e4;
% 
% % 时间戳
% t=0:Ta:N;

f1=5e9;
f2=2e9;

z1=Creat_dither(fs,f1,N);
z2=Creat_dither(fs,f2,N*f2/f1);
