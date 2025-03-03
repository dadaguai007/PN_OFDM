clc;clear;close all;
OFDM_TX;
lw      = 2e6;    % laser linewidth
fs=nn.Fs;
Ta=1/fs;
% fiber param
param=struct();
param.Ltotal = 1000; %km
param.Lspan =10;
param.hz= 10;
param.alpha=0.2;
param.D = 16;
param.gamma = 0;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=fs;

c = 299792458;
c_kms = c/ 1e3 ;   % speed of light (vacuum) in km/s
lamba = c_kms / param.Fc ;

f_deta=fs/(nn.fft_size-1);

ModCarriers=1;
f_k=f_deta*(ModCarriers+nn.nOffsetSub);

f_RF=0;
% 子载波延时
T_k=(param.D *param.Ltotal*lamba.^2)/c_kms*(f_k+f_RF);
% 光载波延时
T_gamma=(param.D *param.Ltotal*lamba.^2)/c_kms*lw;

Ltotal=100:100:5000;
for i=1:length(Ltotal)
    param.Ltotal=Ltotal(i);
    PRT(i)=(2*pi*lw/(Ta*nn.fft_size))*(param.D *param.Ltotal*lamba.^2*f_k/c_kms).^2;
    l(i)=param.Ltotal.^2/nn.fft_size;
end

figure;
plot(Ltotal,PRT)
% figure;
% plot(l)


ModCarriers=1:160;
f_k_squeece=f_deta*(ModCarriers+nn.nOffsetSub);
f_RF=0e9;
% 子载波延时,f_RF越大，则不同子载波的延时改变量越大
T_k_squeece=(param.D *param.Ltotal*lamba.^2)/c_kms*(f_k_squeece+f_RF);
plot(diff(T_k_squeece))
figure;
plot(T_k_squeece)