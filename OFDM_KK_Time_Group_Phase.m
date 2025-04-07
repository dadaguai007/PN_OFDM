% KK 接收 ，分组 Time PN 消除 （失败）

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


% 发射机参数
ofdmPHY=nn;
%%---------------------------------------        解码       ---------------------------%%
Receiver=OFDM_Receiver( ...
    ofdmPHY, ...       %%% 发射机传输的参数
    ofdmPHY.Fs, ...    %   采样
    6*ofdmPHY.Fs, ...  % 上采样
    ofdmPHY.nPkts, ...            % 信道训练长度
    1:4:ofdmPHY.nModCarriers, ...    %导频位置
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
index_carrier=60;
PN_carrier=angle(data_ofdm_martix(index_carrier,:)./qam_signal(index_carrier,:));


% 接收信号中的数据载波
Rec_signal_ofdm=signal_ofdm_martix(Receiver.ofdmPHY.dataCarrierIndex,:);
% 每组数量
Group_Num = 20;
% 每组载波所有载波 进行使用
for m=1:Receiver.ofdmPHY.nModCarriers/Group_Num
    % 每组数据的索引
    Num=(m-1)*Group_Num+1:1:m*Group_Num;
    % 每组所有载波
    Rectr=data_ofdm_martix(Num,:);
    Rec=Rec_signal_ofdm(Num,:);
    H=Hf(Num); % 信道矩阵
    % 信道响应 叠加
    H_data_qam_martix=Rectr.*repmat(H,1,Receiver.ofdmPHY.nPkts);
    % 重建
    ofdm_signal_Rectr=Receiver.Group_Remodulation(H_data_qam_martix,Dc,Group_Num);
    % 接收
    ofdm_signal_Rec=Receiver.Group_Remodulation(Rec,Dc,Group_Num);
    phi_est=angle(ofdm_signal_Rec./ofdm_signal_Rectr);
    % 时域信号
    data=ofdm_signal_Rec.*exp(-1j.*phi_est);
    % 并串转换
    Data(m,:)=data(:);
    Index(m,:)=Num;
end

index=5;
% 转为矩阵形式
ofdm = reshape(Data(index,:),Receiver.ofdmPHY.fft_size,[]);
ofdm = fft(ofdm);
% get the modulated data carriers
data_ofdm = ofdm(Index(index,:),:);
% 信道均衡
% channel estimation
rxTrainSymbol = data_ofdm(:,1:1:Receiver.Nr.nTrainSym);
qam_signal_mat=repmat(Receiver.Implementation.qam_signal,1,1);
refTrainSymbol = qam_signal_mat(Index(index,:),1:1:Receiver.Nr.nTrainSym);
% 信道响应
Hf = mean(rxTrainSymbol./refTrainSymbol,2);

% channel equalization
data_kk = data_ofdm.*repmat(1./Hf,1,Receiver.ofdmPHY.nPkts);
Ref=qam_signal(Index(index,:),:);
% 解码
Receiver.BER(data_kk(:),Ref(:));