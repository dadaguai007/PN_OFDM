classdef OFDM_Receiver < handle

    % 参数直接从 OFDM信号 生成类中 直接引入
    properties
        ofdmPHY; % 传入的信号 参数
        Nr%无量纲参数
        Implementation;% 参考信号实施参数
        Button; % 开关讯号
    end


    methods
        function obj = OFDM_Receiver(varargin)
            if numel(varargin) == 10
                obj.ofdmPHY                     = varargin{1} ;% 传递而来的OFDM 参数
                obj.Nr.fOsc                     = varargin{2};
                obj.Nr.fUp                      = varargin{3};
                obj.Nr.nTrainSym                = varargin{4}; % 训练序列长度
                obj.Nr.pilotIndex               = varargin{5};% 导频位置
                obj.Nr.squ_num                  = varargin{6};% 选取第 x 段信号
                obj.Implementation.ref          = varargin{7}; % 参考序列
                obj.Implementation.qam_signal   = varargin{8};% 调制信号参考矩阵

                obj.Button.CPE_Status           = varargin{9};% 默认 关闭 CPE
                obj.Button.receive_type         = varargin{10};% 默认 直接接收
            else
                error('Number of input variables must be either 0 (default values) or 10');
            end

        end


        function [compSig,lnSig]=KK_receiver(obj,rxSig)
            % 采样倍数
            m=obj.Nr.fUp/obj.Nr.fOsc;
            N=length(rxSig);
            X=fft(rxSig);
            X_up=[X(1:N/2).'  zeros(1,N*m-N)  X(N/2+1:end).'];
            x_up=m*ifft((X_up),m*N);
            upRxSig=x_up;
            %    upRxSig1=resample(rxSig,fUp,fOsc);
            hn = sqrt(upRxSig);
            % ln of abs
            lnSig = log(abs(hn));

            % phi
            H_sig=hilbert(lnSig);

            % exp
            phi = exp(1i.*(imag(H_sig)));
            % phase recovered signal
            compSig = hn.*phi;
        end
        % 信号预处理
        function  ReceivedSignal=Preprocessed_signal(obj,rxsig)
            if strcmp(obj.Button.receive_type,'KK')
                c=0;
                % KK
                [SSB_Sig_KK,ln_sig]=obj.KK_receiver(rxsig+c);

                %下采样
                ReceivedSignal = downsample(SSB_Sig_KK,obj.Nr.fUp/obj.Nr.fOsc);

            elseif strcmp(obj.Button.receive_type,'DD')
                % 不使用KK算法，使用带宽隔开,需要去除DC
                % DC-remove
                rxsig=rxsig-mean(rxsig);
                ReceivedSignal=pnorm(rxsig);
            end

            % 选取某段的信号
            ReceivedSignal=ReceivedSignal(obj.ofdmPHY.len*(obj.Nr.squ_num-1)+1:obj.ofdmPHY.len*obj.Nr.squ_num);
        end


        function [ber,num]=Cal_BER(obj,ReceivedSignal)
            % 解调算法
            [~,~,~,data_qam]=obj.Demodulation(ReceivedSignal);
            % label信号
            ref_seq =qamdemod(ref,obj.ofdmPHY.M,'OutputType','bit','UnitAveragePower',1);
            ref_seq=ref_seq(:);
            % 计算误码率
            [ber,num,~] = CalcBER(data_qam,ref_seq);
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        function [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam]=Demodulation(obj,ReceivedSignal)

            % 解OFDM
            signal_ofdm = reshape(ReceivedSignal, obj.ofdmPHY.fft_size+ obj.ofdmPHY.nCP,[]); % 转换为矩阵形式
            signal_ofdm(1: obj.ofdmPHY.nCP,:) = [];% 去除CP
            signal_ofdm = fft(signal_ofdm);
            % 存储矩阵形式的OFDM 信号
            signal_ofdm_martix=signal_ofdm;
            % get the modulated data carriers
            data_ofdm = signal_ofdm(obj.ofdmPHY.postiveCarrierIndex,:);
            % 信道均衡
            [data_ofdm,Hf]=obj.one_tap_equalization(data_ofdm);

            % CPE compensation
            if strcmp(obj.Button.CPE_Status,'on')
                data_ofdm=CPE_Eliminate(data_ofdm);
            end

            %保留信号矩阵
            data_ofdm_martix=data_ofdm;
            %归一化
            data_ofdm=data_ofdm(:);% 矩阵转换为行向量
            data_ofdm = data_ofdm./sqrt(mean(abs(data_ofdm(:)).^2));
            % 硬判决
            data_qam=hard_decision(data_ofdm);

        end

        function ofdmSig=Remodulation(obj,ReceivedSignal)
            % 解调获得 信号
            [~,~,~,data_qam]=obj.Demodulation(ReceivedSignal);
            % 重新调制
            % nOffsetSub 行置零 ,positive 放置qam ， 后续置零
            X= ([zeros(obj.ofdmPHY.nOffsetSub,obj.ofdmPHY.nPkts);...
                data_qam; ...
                zeros(obj.ofdmPHY.fft_size-obj.ofdmPHY.nModCarriers-obj.ofdmPHY.nOffsetSub,obj.ofdmPHY.nPkts)]);
            % 转换为时域
            ofdmSig=ifft(X);
            % 添加CP
            ofdmSig = [ofdmSig(end-obj.ofdmPHY.nCP+1:end,:);ofdmSig];
            % 并串转换
            ofdmSig = ofdmSig(:);
        end


        % ZF信道均衡
        function [data_kk,Hf]=one_tap_equalization(obj,data_ofdm)
            % channel estimation
            rxTrainSymbol = data_ofdm(:,1:obj.Nr.nTrainSym);
            qam_signal_mat=repmat(obj.Implementation.qam_signal,1,1);
            refTrainSymbol = qam_signal_mat(:,1:obj.Nr.nTrainSym);
            % 信道响应
            Hf = mean(rxTrainSymbol./refTrainSymbol,2);

            % channel equalization
            data_kk = data_ofdm.*repmat(1./Hf,1,obj.Nr.nPkts);

        end
        % 相位均衡
        function data=CPE_Eliminate(obj,data_ofdm)

            phi_mean=angle(mean(data_ofdm(obj.Nr.pilotIndex,:)./...
                obj.Implementation.qam_signal(obj.Nr.pilotIndex,:),1));


            data=data_ofdm.*...
                repmat(exp(-1j.*phi_mean),size(data_ofdm,1),1);
        end
        % 硬判决
        function data_qam=hard_decision(obj,Receiver_data)
            data_qam = qamdemod(Receiver_data,obj.ofdmPHY.M,'OutputType','bit','UnitAveragePower',1);
            data_qam=data_qam(:);
        end
        % EVM测算

        function [rmsEVM_symbol,rmsEVM_subcarrier,rmsEVM_martix]=EVM_Measure_martix(obj,data_ofdm_martix)
            %EVM measure
            evm = comm.EVM(AveragingDimensions=1);
            rmsEVM_symbol = evm(data_ofdm_martix,obj.Implementation.qam_signal);


            evm = comm.EVM(AveragingDimensions=2);
            rmsEVM_subcarrier = evm(data_ofdm_martix,obj.Implementation.qam_signal);


            evm = comm.EVM(AveragingDimensions=[1 2]);
            rmsEVM_martix = evm(data_ofdm_martix,obj.Implementation.qam_signal);

            fprintf('EVM = %1.7f\n',rmsEVM_martix);
        end

        function rmsEVM_subcarrier=EVM_Measure_symbol(obj,data_ofdm_martix,symbol_indedx)
            %EVM measure
            % 某个符号上所有载波的EVM
            evm = comm.EVM(AveragingDimensions=2);
            rmsEVM_subcarrier = evm(data_ofdm_martix(:,symbol_indedx),obj.Implementation.qam_signal(:,symbol_indedx));

        end

        % 星座图绘制
        function scatter_plot(obj,ReceivedSignal)
            [~,data_ofdm_martix,~,~]=obj.Demodulation(ReceivedSignal);
            dara_ofdm_squ=data_ofdm_martix(:);
            scatterplot(dara_ofdm_squ);
        end

    end


end