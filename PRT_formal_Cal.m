clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')
% 使用公式进行计算：参考文章Analysis of laser phase noise effect in direct- detection optical OFDM transmission

for index_f=1:5

OFDM_TX_subcarrier_Add;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);

% System Parameters
fs=nn.Fs;
Ta=1/fs;
%
scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal_ofdm = signal./scale_factor;

% Carrier
A=1;
%Amp_factor
Amp=1;


% signal_TX
signal_TX=A+Amp*signal_ofdm.';


% 记录最原始的数据长度
N_index=length(signal_TX);

%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);

% Laser
% Generate LO field with phase noise
% 输入光功率
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
Ai= sqrt(Pi);
lw      = 2e6;    % laser linewidth
phi_pn_lo = phaseNoise(lw, length(signal_TX), Ta);
sigLO = exp(1j * phi_pn_lo);
Pin=Ai*sigLO;

% Optical Signal
signal_TXO=signal_TX.*Pin;
% signal_TXO=signal_TX;
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

% Power_PRT Cal
c = 299792458;
c_kms = c/ 1e3 ;   % speed of light (vacuum) in km/s
lamba = c_kms / param.Fc ;
T=f_idex*lamba.^2*param.Ltotal*param.D/c_kms;
M_k=T/Ta;
PRT_power(index_f)=(2*pi*lw*Ta*M_k/nn.fft_size.^2)*(-1/3*M_k.^2+M_k*nn.fft_size+1/3);


%  ICI power Cal
beta_k=2*pi*lw*T;
ICI_up_power=(beta_k/nn.fft_size.^2)*(nn.fft_size.^2-nn.fft_size*M_k+1/3*M_k.^2-1/3);


end

figure;
plot(PRT_power)