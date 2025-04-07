clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
% 发射机配置(需要配置为KK模式)
OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);


% 归一化
scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal_ofdm = signal./scale_factor;


% 转置
signal_ofdm=signal_ofdm.';


% System Parameters
% 上采样
m=5;
fs=m*nn.Fs;
Ta=1/fs;
% 信号上采样
signal_ofdm=nn.Up_sample(m,signal_ofdm.');


%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);

%%---------------------------------------          Laser           ----------------------------%%
% 时间轴
[~,t_up]=freq_time_set(length(signal_ofdm),fs);
% Generate LO field with phase noise
% 输入光功率
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
lw      = 10e6;    % laser linewidth
phi_pn_lo = phaseNoise(lw, length(signal_ofdm), Ta);
sigLO = exp(1j * phi_pn_lo);
% Lasear 1
f_lo1      = 0;                 % frequency offset
sigLO1 = sigLO.*exp(1j*2*pi*f_lo1*t_up); % add frequency offset
% Lasear 2
f_lo2      = 100e9;                 % frequency offset
sigLO2 = sigLO.*exp(1j*2*pi*f_lo2*t_up); % add frequency offset

%Amp_factor
Amp=1.8;

%%---------------------------------------          Modulator           ----------------------------%%
% Optical Signal
nn.ModulationPHY.Pi_dBm=Pi_dBm; %dB
nn.ModulationPHY.Amp=Amp; % 信号放大
nn.ModulationPHY.Vpi=10; % Vpi
bias_phi=0.89;
% 调制
signal_TXO1=nn.OFDM_Modulation(bias_phi,signal_ofdm,sigLO1);
% CSPR
CSPR1=nn.Cal_CSPR(signal_TXO1,bias_phi,sigLO1);


% Optical Signal
nn.ModulationPHY.Pi_dBm=Pi_dBm; %dB
nn.ModulationPHY.Amp=Amp; % 信号放大
nn.ModulationPHY.Vpi=10; % Vpi
bias_phi=0.89;
% 调制
signal_TXO2=nn.OFDM_Modulation(bias_phi,signal_ofdm,sigLO2);
% CSPR
CSPR2=nn.Cal_CSPR(signal_TXO2,bias_phi,sigLO2);

% 信号合并
signal_TXO=signal_TXO1+signal_TXO2;


% fiber param
param=struct();
param.Ltotal = 80; %km
param.Lspan =10;
param.hz= 1;
param.alpha=0.2;
param.D = 16;
param.gamma = 0;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=fs;



%noise
%sigTxo=awgn(sigTxo,snr(index),'measured');

% Transmission
Train_type='ssfm';
if strcmp(Train_type,'ssfm')
    OFDM_Train;
else
    sigRxo=signal_TXO;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           接收         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% 采样率恢复正常

% 光滤波器 1
dataout1=LPF(sigRxo,fs,fs/m/2);
sigRxo1=downsample(dataout1,m);


% 光滤波器 2
sigRxo_Wave_Select=sigRxo.*exp(-1j*2*pi*f_lo2*t_up);
dataout2=LPF(sigRxo_Wave_Select,fs,fs/m/2);
sigRxo2=downsample(dataout2,m);


% PD param
paramPD=struct();
paramPD.B =fs/m;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs/m;
% DD-OFDM
%pd-Receiver
ipd_btb = pd(sigRxo1, paramPD);
mon_ESA(ipd_btb,fs/m);

% 发射机参数
ofdmPHY=nn;
%%---------------------------------------        解码       ---------------------------%%
Receiver=OFDM_Receiver( ...
                        ofdmPHY, ...       %%% 发射机传输的参数
                        ofdmPHY.Fs, ...    %   采样
                        6*ofdmPHY.Fs, ...  % 上采样
                        ofdmPHY.nPkts, ...            % 信道训练长度 1:1:ofdmPHY.nPkts
                        1:1:ofdmPHY.nModCarriers, ...    %导频位置
                        1, ...             % 选取第一段信号
                        ref_seq, ...       % 参考序列
                        qam_signal, ...    % qam 矩阵
                        'off', ...         % 是否采用CPE
                        'off', ...         % 对所有载波进行相位补偿
                        'KK');             % 接收方式

% 信号预处理
[ReceivedSignal,~]=Receiver.Preprocessed_signal(ipd_btb);
% 归一化
ReceivedSignal=pnorm(ReceivedSignal);
% 直流
Dc=mean(ReceivedSignal);
% BER 计算
[ber,num]=Receiver.Cal_BER(ReceivedSignal);

%  解调
[signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Receiver.Demodulation(ReceivedSignal);

% 相噪
index_carrier=110;
PN_carrier=angle(data_ofdm_martix(index_carrier,:)./qam_signal(index_carrier,:));
