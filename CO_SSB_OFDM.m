clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);

% System Parameters
fs=nn.Fs;
Ta=1/fs;
% 归一化
scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal_ofdm = signal./scale_factor;

% Carrier
A=1;
%Amp_factor
Amp=1;
% 基带OFDM信号转置
signal_ofdm=signal_ofdm.';


% signal_TX
signal_TX=A+Amp*signal_ofdm;


% CSPR measure
CSPR_Mea;

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
lw      = 100e3;    % laser linewidth
phi_pn_lo = phaseNoise(lw, length(signal_TX), Ta);
sigLO = exp(1j * phi_pn_lo);
Pin=Ai*sigLO;


% Optical Signal
signal_TXO=signal_TX;


% fiber param
param=struct();
param.Ltotal = 700; %km
param.Lspan =10;
param.hz= 1;
param.alpha=0.2;
param.D = 16;
param.gamma = 0;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=fs;

% PD param
paramPD=struct();
paramPD.B =fs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs;



% 重复信号
k=1;
% qam信号矩阵
ref_seq_mat=repmat(qam_signal,1,k);


%noise
%sigTxo=awgn(sigTxo,snr(index),'measured');

%power
power=signalpower(signal_TXO);
fprintf('optical signal power: %.2f dBm\n', 10 * log10(power / 1e-3));


% Transmission
Train_type='ssfm';
if strcmp(Train_type,'ssfm')
    OFDM_Train;
else
    sigRxo=signal_TXO;
end


% CO-OFDM
% linear-Receiver
% phase noise
ipd_btb=sigRxo.*sigLO;

% CPE compensation 参数
CPE_Status='off';
if strcmp(CPE_Status,'on')
    W=160;
    pilotIndex=1:2:W;
end


% 按顺序解码
OFDM_Decode_squence;

% Rx_16QAM
data_kk_mat;


% 第i个载波上的相噪分布
jj=1;
phase_jj=angle(data_kk_mat(jj,:)./qam_signal_mat(jj,:));

% PRT方差
phase_var=var(phase_jj);
% 功率
PRT_power=signalpower(phase_jj);
fprintf('PRT_power = %1.7f\n',PRT_power);

figure;
plot(phase_jj)
xlabel('symbol')
ylabel('phas fluctuation/rad')


scatterplot(signal_scatter);