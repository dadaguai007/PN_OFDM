% PCP 方式 ， 采用KK 探测
clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

% 信号发射机 KK  配置
OFDM_TX_phase_conjugated;
% 生成信号
[y1,y2,signal,qam_signal,index_data,index_pcp,qam_mat_phase]=nn.Output();
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

% 转置
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
lw      = 10e6;    % laser linewidth
phi_pn_lo = phaseNoise(lw, length(signal_TX), Ta);
sigLO = exp(1j * phi_pn_lo);
Pin=Ai*sigLO;


% Optical Signal
signal_TXO=signal_TX.*sigLO;

% fiber param
param=struct();
param.Ltotal = 700; %km
param.Lspan =10;
param.hz= 1;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
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

% DD-OFDM
%pd-Receiver
ipd_btb = pd(sigRxo, paramPD);
mon_ESA(ipd_btb,fs);


% 发射机参数
ofdmPHY=nn;
%%---------------------------------------        解码       ---------------------------%%
Receiver=OFDM_Receiver( ...
                        ofdmPHY, ...       %%% 发射机传输的参数
                        ofdmPHY.Fs, ...    %   采样
                        6*ofdmPHY.Fs, ...  % 上采样
                        ofdmPHY.nPkts, ...            % 信道训练长度
                        1:1:ofdmPHY.nModCarriers, ...    %导频位置
                        1, ...             % 选取第一段信号
                        ref_seq, ...       % 参考序列
                        qam_mat_phase, ...    % qam_pcp 矩阵
                        'off', ...         % 是否采用CPE
                        'off', ...         % 对所有载波进行相位补偿
                        'KK');             % 接收方式

% 信号预处理
[ReceivedSignal,Dc]=Receiver.Preprocessed_signal(ipd_btb);
% BER 计算
[ber,num]=Receiver.PCP_Cal_BER(ReceivedSignal,index_data,index_pcp);

