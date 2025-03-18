% OFDM_Decode_squence PCP  qam信号与 qam——共轭各一半
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
    qam_signal_mat=repmat(qam_mat_phase,1,HK);
    refTrainSymbol = qam_signal_mat(:,1:nTrainSym);
    Hf = mean(rxTrainSymbol./refTrainSymbol,2);

    % channel equalization
    data_kk = data_kk_ofdm.*repmat(1./Hf,1,nn.nPkts*HK);

    % CPE compensation
    if strcmp(CPE_Status,'on')
        phase_compensation;
    end

    % 相位共轭进行处理
    % 找到索引
    data_index=index_data-nn.nOffsetSub;
    pcp_index=index_pcp-nn.nOffsetSub;
    % 选取共轭数据
    data_kk_data=data_kk(data_index,:);
    data_kk_pcp=data_kk(pcp_index,:);


    %   共轭消除

    data_rec=(data_kk_data+conj(data_kk_pcp))/2;
    %保留信号矩阵
    data_kk_mat=data_rec;
    %归一化
    data_rec=data_rec(:);
    data_rec = data_rec./sqrt(mean(abs(data_rec(:)).^2));


    % 解码，计算误码
    % 参考序列
    ref_seq_1 =qamdemod(ref_seq,nn.M,'OutputType','bit','UnitAveragePower',1);
    ref_seq_1=ref_seq_1(:);
    % 接收序列
    yyy = data_rec;

    yyy_1 = qamdemod(yyy,nn.M,'OutputType','bit','UnitAveragePower',1);
    yyy_1=yyy_1(:);
    % 全部的序列进行解码
    [ber1(squ_num),num1(squ_num),error_location] = CalcBER(yyy_1,ref_seq_1); %计算误码率
    fprintf('Num of Errors = %d, BER = %1.7f\n',num1(squ_num),ber1(squ_num));
    % 按符号解码，每个符号上的错码数
    Calc_BER_mat;


end
