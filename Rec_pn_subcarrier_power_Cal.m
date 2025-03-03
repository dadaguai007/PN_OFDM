clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')
% 仿真分析Total phase noise 、PRT 、ICI的分布
%% PRT

% 设置NUM个载波进行PRT和ICI计算
NUM=50;
% 确定信号长度
OFDM_TX;
% 生成信号
[~,~,signal,~,~]=nn.Output();
% System Parameters
fs=nn.Fs;
Ta=1/fs;
% 相噪固定
% Laser
% Generate LO field with phase noise
% 输入光功率
Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W
Ai= sqrt(Pi)*10;
lw      = 4e6;    % laser linewidth
phi_pn_lo = phaseNoise(lw, length(signal), Ta);
sigLO = exp(1j * phi_pn_lo);
Pin=Ai*sigLO;



% NUM=30;
PTR=zeros(NUM,1);
for index_f=1:NUM

    OFDM_TX_subcarrier_Add;
    % 生成信号
    [~,~,signal,qam_signal,postiveCarrierIndex]=nn.Output();
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



    % Optical Signal

    signal_TXO=signal_TX.*sigLO;


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

    % DD-OFDM
    %pd-Receiver
    ipd_btb = pd(sigRxo, paramPD);
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
    fprintf('Recceivver_power = %1.7f\n',phase_var);

    PTR(index_f)=phase_var;


end
%%
% ICI Cal and Total Phase noise

% 载波数据填满传输情况，用来计算ICI
OFDM_TX_single_subcarrier;
% 生成信号
[~,~,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);

% System Parameters
fs=nn.Fs;
Ta=1/fs;
%
scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal_ofdm = signal./scale_factor;

signal_ofdm=signal_ofdm.';



% signal_TX
signal_TX=A+Amp*signal_ofdm;


% 记录最原始的数据长度
N_index=length(signal_TX);

%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);


% Optical Signal
signal_TXO=signal_TX.*sigLO;


% 重复信号
k=1;
% qam信号矩阵
ref_seq_mat=repmat(qam_signal,1,k);

%power
power=signalpower(signal_TXO);
fprintf('optical signal power: %.2f dBm\n', 10 * log10(power / 1e-3));


% Transmission
OFDM_Train;


% DD-OFDM
%pd-Receiver
ipd_btb = pd(sigRxo, paramPD);

% CPE compensation 参数
CPE_Status='off';
if strcmp(CPE_Status,'on')
    W=160;
    pilotIndex=1:2:W;
end

% 按顺序解码
OFDM_Decode_matrix;

% Rx_16QAM
data_kk_mat;


% 第i个载波上的总相噪分布
for index=1:index_f
    phase_total(index,:)=angle(data_kk_mat(index,:)./qam_signal_mat(index,:));
    % total phase power
    phase_total_var(index)=var(phase_total(index,:));
end
% Cal ICI
ICI=phase_total_var.'-PTR;

% ICI的直方图
histogram(ICI)
% PTR 曲线
figure;hold on;
plot(PTR)
plot(ICI)
xlabel('subcarrier')
ylabel('power')