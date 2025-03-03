clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')

% 复现：Laser phase noise and OFDM symbol duration effects on the performance of direct-detection based optical OFDM access network

% 不同长度下，单载波传输和多载波传输的PN相噪大小
L=200:200:5000;

PTR=zeros(length(L),1);

% 假设传输Num个载波
% 计算第index_f个载波的PRT、ICI、PN
index_f=1;
Num=30;
for index=1:length(L)
    % 单载波传输
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

    % 相噪固定
    % Laser
    % Generate LO field with phase noise
    % 输入光功率
    Pi_dBm = 10;
    Pi = 10^(Pi_dBm/10)*1e-3; %W
    Ai= sqrt(Pi)*10;
    lw      = 5e6;    % laser linewidth
    phi_pn_lo = phaseNoise(lw, length(signal), Ta);
    sigLO = exp(1j * phi_pn_lo);
    Pin=Ai*sigLO;

    % Optical Signal

    signal_TXO=signal_TX.*sigLO;


    % fiber param
    param=struct();
    param.Ltotal = L(index); %km
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

    PTR(index)=phase_var;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 所有数据载波的传输，计算ICI和PN_total

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

    % 第i个载波上的总相噪分布,对应单载波信号传输的载波位置

    phase_total=angle(data_kk_mat(jj,:)./qam_signal_mat(jj,:));
    % total phase power
    phase_total_var(index)=var(phase_total);
    % Cal ICI
    ICI(index)=phase_total_var(index)-PTR(index);


end



% PTR 曲线
figure;hold on;
plot(L,PTR)
% plot(L,ICI)
% plot(L,phase_total_var)
xlabel('Transmission length(km)')
ylabel('power')
legend('PRT','ICI','PN')