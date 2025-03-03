clc;clear;close all;
% 复现：Analytical study of optical SSB-DMT with IMDD
addpath('Fncs\')
addpath('D:\PhD\Codebase\')

% OFDM 信号
OFDM_TX;

fs=nn.Fs;
Ta=1/fs;
% 频率间隔
f_deta=fs/(nn.fft_size-1);

% fiber param
param=struct();
param.Ltotal = 300; %km
param.Lspan =10;
param.hz= 10;
param.alpha=0.2;
param.D = 16;
param.gamma = 0;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=fs;

lw      = 100e3;    % laser linewidth

c = 299792458;
c_kms = c/ 1e3 ;   % speed of light (vacuum) in km/s
lamba = c_kms / param.Fc ; % waveform


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P2A Cal （根据公式进行计算）
% syms f
k_dot=lamba.^2*param.D/(2*pi*c_kms);
% besselTerm=@(f) 4*besselj(n,1/f*sqrt(2*lw/pi))*besselj(n+1,1/f*sqrt(2*lw/pi));
% sineTerm=@(f) sin(1/2*(2*n+1)*(2*pi*f).^2*k_dot*param.Ltotal);
% Term=@(f) 4*besselj(n,1/f*sqrt(2*lw/pi))*besselj(n+1,1/f*sqrt(2*lw/pi))*...
%     sin(1/2*(2*n+1)*(2*pi*f).^2*k_dot*param.Ltotal);

num=1:512;
for j=1:length(num)

    % 载波数
    index=num(j);
    % index=50;
    F_deta=linspace((index-1)*f_deta,index*f_deta,100);

    for i=1:length(F_deta)
        f=F_deta(i);
        n=0:1:100;
        besselTerm= 4*besselj(n,1/f*sqrt(2*lw/pi)).*besselj(n+1,1/f*sqrt(2*lw/pi));
        sineTerm= sin(1/2*(2*n+1)*((2*pi*f).^2)*k_dot*param.Ltotal);
        x(i)=(1/2)*(sum(besselTerm.*sineTerm).^2);
    end

    power(j)=sum(x);
end

% 信号载波功率比
gamma=[-30.5;-27;-22.5];
T=10.^(gamma/10);

ratioTerm=repmat(T,1,length(power))./nn.nModCarriers;

% P2A 的SNR
SNR_P2A=ratioTerm./power;
SNR_P2A_dB=10*log10(SNR_P2A);
figure;hold on;box on;
plot(num,SNR_P2A_dB(1,:))
plot(num,SNR_P2A_dB(2,:))
plot(num,SNR_P2A_dB(3,:))
title('SNR-P2A')
xlabel('subcarrier')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSBI（无法计算）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICI和PRT计算（根据公式计算）

k=1:512;
% 第k子信道的时间延迟或离场
T_k=k*c_kms*param.D*param.Ltotal*f_deta/(param.Fc.^2);

beta_k=2*pi*lw*T_k;
%第k个子信道相对于载波的延迟（采样数）
M_k=T_k/Ta;
% 功率衰减
alpha_k=1-beta_k;

% ICI power
coefficientTerm=(beta_k/(nn.fft_size.^2));
factorTerm=M_k.^2/3+nn.fft_size.^2-M_k*nn.fft_size-1/3;
ICI_power=coefficientTerm.*factorTerm;

% PRT power
PRcoefficientTerm=(beta_k/(3*(nn.fft_size.^2)));
PRfactorTerm=-M_k.^2+3*M_k*nn.fft_size+1;
PRT_power=PRcoefficientTerm.*PRfactorTerm;

% ICI SNR

SNR_ICI_dB=10*log10(alpha_k./ICI_power);

figure;
plot(SNR_ICI_dB)
title('SNR-ICI')
xlabel('subcarrier')

figure;hold on;
plot(k,PRT_power)
plot(k,ICI_power)
legend('PRT','ICI')
xlabel('subcarrier')