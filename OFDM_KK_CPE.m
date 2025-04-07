% KK 接收 ，CPE 和 total 载波共同消除

clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
% 发射机配置(需要配置为KK模式)
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


% 转置
signal_ofdm=signal_ofdm.';


%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);

%%---------------------------------------          Laser           ----------------------------%%
% Generate LO field with phase noise
% 输入光功率
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
Ai= sqrt(Pi);
lw      = 10e6;    % laser linewidth
phi_pn_lo = phaseNoise(lw, length(signal_ofdm), Ta);
sigLO = exp(1j * phi_pn_lo);


%Amp_factor
Amp=1.8;

%%---------------------------------------          Modulator           ----------------------------%%
% Optical Signal
nn.ModulationPHY.Pi_dBm=Pi_dBm; %dB
nn.ModulationPHY.Amp=Amp; % 信号放大
nn.ModulationPHY.Vpi=10; % Vpi
phi=0.89;
% 调制
signal_TXO=nn.OFDM_Modulation(phi,signal_ofdm,sigLO);
% CSPR
CSPR=nn.Cal_CSPR(signal_TXO,phi,sigLO);
% fiber param
param=struct();
param.Ltotal = 200; %km
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

% 采用CPE补偿方式 测试

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
                        qam_signal, ...    % qam 矩阵
                        'on', ...         % 是否采用CPE
                        'off', ...         % 对所有载波进行相位补偿
                        'KK');             % 接收方式

% 信号预处理
[ReceivedSignal,Dc]=Receiver.Preprocessed_signal(ipd_btb);
% BER 计算
[ber,num]=Receiver.Cal_BER(ReceivedSignal);

%  解调
[signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Receiver.Demodulation(ReceivedSignal);

% 相噪
index_carrier=60;
PN_carrier=angle(data_ofdm_martix(index_carrier,:)./qam_signal(index_carrier,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 对所有载波进行相位消除       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Receiver.Button.CPE_Status           = 'off';% 默认 关闭 CPE
Receiver.Button.PN_Total_Carrier     = 'on';% 打开 所有载波相除相噪
% BER 计算
[ber1,num1]=Receiver.Cal_BER(ReceivedSignal);

figure;
plot(PN_carrier)
title('第60个载波的相噪分布')