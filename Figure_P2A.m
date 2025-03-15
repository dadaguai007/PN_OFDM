clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

 
% % System Parameters
fs=40e9;
Ta=1/fs;


% 数据长度
data_num=1e6;
% Laser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Generate LO field with phase noise %%%%%%%%%%%%%%%%%%%%%%%%
% 输入光功率
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
Ai= sqrt(Pi);
lw      = 2e6;    % laser linewidth
phi_pn_lo = phaseNoise(lw, data_num, Ta);
sigLO = exp(1j * phi_pn_lo);
Pin=Ai*sigLO;


% fiber param
param=struct();
param.Ltotal = 80; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='none';
param.Fs=fs;

% PD param
paramPD=struct();
paramPD.B =fs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs;




%power
power=signalpower(Pin);
fprintf('optical signal power: %.2f dBm\n', 10 * log10(power / 1e-3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           P2A Train         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigRxo=ssfm(Pin,param);
power2=signalpower(sigRxo);
fprintf(' after ssfm signal power: %.2f dBm\n', 10 * log10(power2 / 1e-3));

%pd-Receiver
ipd_btb = pd(sigRxo, paramPD);
mon_ESA(ipd_btb,fs);


