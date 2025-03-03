% OFDM_Decode_matrix
% 只对OFDM信号进行解码，整体矩阵进行EVM测算，不进行按照符号之类的性能测算
rxsig = ipd_btb(1:2*floor(length(ipd_btb)/2)).';

receive_type='DD';
if strcmp(receive_type,'KK')
    c=0;
    % KK
    fup=fs*4;
    [SSB_Sig_KK,ln_sig] = KK_MPC_Multi(rxsig+c,fs,fup);

    %下采样
    data_kk = downsample(SSB_Sig_KK,fup/fs);
    DATA=data_kk;
elseif strcmp(receive_type,'DD')
    % 不使用KK算法，使用带宽隔开,需要去除DC
    % DC-remove
    rxsig=rxsig-mean(rxsig);
    DATA=pnorm(rxsig);
end


for squ_num=1:k
    % 训练序列
    nTrainSym=10;
    HK=1;
    % 选取每一个阶段的信号
    DATA_squ=DATA(N_index*(squ_num-1)+1:N_index*squ_num);


    % 解OFDM 信道均衡
    data_kk_ofdm = reshape(DATA_squ,nn.fft_size+nn.nCP,[]);
    data_kk_ofdm(1:nn.nCP,:) = [];
    data_kk_ofdm = fft(data_kk_ofdm);
    Data_matrix=data_kk_ofdm;% 存储RX_OFDM矩阵
    % get the modulated carriers
    data_kk_ofdm = data_kk_ofdm(postiveCarrierIndex,:);
    % channel estimation
    rxTrainSymbol = data_kk_ofdm(:,1:nTrainSym);
    qam_signal_mat=repmat(qam_signal,1,HK);
    refTrainSymbol = qam_signal_mat(:,1:nTrainSym);
    Hf = mean(rxTrainSymbol./refTrainSymbol,2);

    % channel equalization
    data_kk = data_kk_ofdm.*repmat(1./Hf,1,nn.nPkts*HK);

    % CPE compensation
    if strcmp(CPE_Status,'on')
        phase_compensation;
    end

    
    %保留信号矩阵
    data_kk_mat=data_kk;
    %归一化
    data_kk=data_kk(:);
    data_kk = data_kk./sqrt(mean(abs(data_kk(:)).^2));

    % 解码，计算误码
    % 参考序列
    ref_seq_1 =qamdemod(ref_seq,nn.M,'OutputType','bit','UnitAveragePower',1);
    ref_seq_1=ref_seq_1(:);
    % 接收序列
    yyy = data_kk;

    yyy_1 = qamdemod(yyy,nn.M,'OutputType','bit','UnitAveragePower',1);
    yyy_1=yyy_1(:);
    % 全部的序列进行解码
    [ber1(squ_num),num1(squ_num),error_location] = CalcBER(yyy_1,ref_seq_1); %计算误码率
    fprintf('Num of Errors = %d, BER = %1.7f\n',num1(squ_num),ber1(squ_num));
    % 按符号解码，每个符号上的错码数
    Calc_BER_mat;

    OFDM_EVM;
%     symbol_EVM(:,squ_num)=rmsEVM_symbol.';
%     subcarrier_EVM(:,squ_num)=rmsEVM_subcarrier;
% 
%     % 第L个符号的EVM
%     LL=200;
%     OFDM_symbol_EVM;
%     subcarrier_index_symbol_EVM(:,squ_num)=rmsEVM_subcarrier_index_symbol;
% 
%     % Rx_16QAM
%     signal_scatter(:,squ_num)=data_kk;
%     %Rx_OFDM
%     signal_squ_ofdm(:,squ_num)=DATA_squ;

end
