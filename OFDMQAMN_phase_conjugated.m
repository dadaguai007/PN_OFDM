classdef OFDMQAMN_phase_conjugated < handle
    % OFDM信号的生成 ，包括 脉冲成型步骤 ，以及频率谱的绘制
    % 输出信号：  SSB-PAM 信号
    % 随机码包括：PRBS ，rand 两种
    properties
        order;
        prbsOrder;
        M;
        NSym; %码元数目
        DataType ;
        Rs; % baud rate
        Fs; % sample rate
        Nsam;%每个码元周期的采样次数

        psfShape;
        psfshape;
        psfRollOff;
        psfLength;

        ipaddr
        signalname
        defaul

        fft_size;
        nPkts;
        nCP;
        nModCarriers;
        nOffsetSub;

        percentage; % PCP的百分比
        Type

        len;
        dataCarrierIndex;
        ModulationPHY;%信号调制参数
    end

    properties (Dependent = true)
        hsqrt;
        psf;

    end

    methods
        function obj=OFDMQAMN(IP)
            if nargin < 1
                obj.ipaddr = 'no value';
            else
                obj.ipaddr = IP;
            end
            obj.M =4;
            obj.Rs=8e9;
            obj.NSym = 32768;
            obj.Fs = 64e9;
            obj.prbsOrder = 15;
            obj.Nsam = obj.Fs/obj.Rs;
            obj.psfShape = 'sqrt';
            obj.psfRollOff = 0.01;
            obj.psfLength=256;
            obj.fft_size = 1024;
            obj.nPkts = 100;
            obj.nCP = 32;
            obj.nModCarriers = 320;
            obj.percentage=1/3; % 默认
            obj.Type ='std';
        end

        %生成表头文件
        function summaryTable = summary(obj)
            disp('Nyquist PAM-N class parameters summary:');
            rows = {'PRBS order';'PAM order';'Number of symbols';'Baud rate';...
                'Sampling rate';'Pulse shaping filter shape';'Roll off';...
                'Length of filter';'Sample per symbol'};
            varName = {'prbsOrder';'M';'NSym';'Rs';'Fs';'psfShape';...
                'psfRollOff';'psfLength';'Nsam'};
            values = {obj.prbsOrder;obj.M;obj.NSym;obj.Rs/1e9;obj.Fs/1e9;obj.psfShape;...
                obj.psfRollOff;obj.psfLength;obj.Nsam};
            units = {' ';' ';' ';'Gbaud';'GSample/s';' ';' ';' ';' '};

            summaryTable = table(varName,values,units,'RowNames',rows);
            disp(summaryTable);
        end

        % 生成 索引  与 qam组合矩阵
        function [y1,y2,signal,qam_signal,index_data,index_pcp,qam_mat_phase]=Output(obj)
            %生成bit数，进行pammod的调制,脉冲成型，希尔伯特变换，取实部和虚部一起输出
            switch lower(obj.DataType)
                case 'prbs'
                    [~,symbols]=obj.prbs_bits();
                case 'rand'
                    symbols=obj.rand_bits();
            end
            %调制
            qam_signal=obj.qam(symbols);
            %归一化
            %qam_signal = qam_signal./sqrt(mean(qam_signal.^2));

            % 50%PCP
            [ofdm_signal_pcp,index_data,index_pcp,qam_mat_phase] = obj.ofdm_pcp_50(qam_signal);

            %                [ofdm_signal_pcp,index_data,index_pcp] = obj.ofdm_pcp_33(qam_signal);
            %重采样
            signal =resample(ofdm_signal_pcp,obj.Fs,obj.Rs);
            SigI=real(signal);
            SigQ=imag(signal);
            y1=SigI(:);
            y2=SigQ(:);
        end
% 信号调制
        function sigTxo=OFDM_Modulation(obj,phi,signal,Lo)
            % IQ调制
            % 输入光功率
            % 转换为W
            Pi=10^(obj.ModulationPHY.Pi_dBm/10)*1e-3; %W
            % phi 为 Vbias 偏移程度   phi=0.87; 
            Amp=obj.ModulationPHY.Amp; % 信号放大      Amp=1.7;
            % 调制器参数
            param.Vpi=obj.ModulationPHY.Vpi;
            param.VbI=-phi*obj.ModulationPHY.Vpi;
            param.VbQ=-phi*obj.ModulationPHY.Vpi;
            param.Vphi=param.Vpi/2;
            % 输入光功率
            Ai  = sqrt(Pi).*Lo;
            % 调制指数
            m=Modulation_index(Amp*signal.',param.Vpi,'ofdm');
            fprintf('Modulation index=%3.3f\n',m);
            % 调制
            sigTxo = iqm(Ai, Amp*signal, param);
            % 功率计算
            signal_power=signalpower(sigTxo);
            fprintf('optical signal power: %.2f dBm\n', 10 * log10(signal_power / 1e-3));
        end
        % CSPR 计算
        function CSPR=Cal_CSPR(obj,sigTxo,phi,Lo)

            % Cal CSPR
            param.Vpi=obj.ModulationPHY.Vpi;
            param.VbI=-phi*obj.ModulationPHY.Vpi;
            param.VbQ=-phi*obj.ModulationPHY.Vpi;
            param.Vphi=obj.ModulationPHY.Vpi/2;
             % 输入光功率
             Pi=10^(obj.ModulationPHY.Pi_dBm/10)*1e-3; %W
            Ai  = sqrt(Pi);
            Eout1 = iqm(Ai*Lo, 0, param);

            Eout_s=sigTxo-Eout1;
            power1=signalpower(Eout1);
            Ps=signalpower(Eout_s);
            CSPR = power1/(Ps);
            CSPR=10*log10(CSPR);
            fprintf("the CSPR is %1.2fdB\n",CSPR);

        end
        function [data,symbols_prbs]=prbs_bits(obj)
            %参数：obj.prbsOrder，NSym，M
            %采用prbs码生成基本数据
            data = prbs1(obj.prbsOrder,obj.NSym*log2(obj.M),0);
            data_2bit=reshape(data,log2(obj.M),[]);
            %                         symbols_prbs1 = 2.^(0:log2(obj.M)-1)*data_2bit;
            %             symbols_prbs = reshape(symbols_prbs,obj.nModCarriers,[]);
            symbols_prbs=double(data_2bit);

            %进行重组，映射
            %             data_2bit=reshape(data,length(data)/log2(obj.M),[]);
            %             data_2bit = data_2bit.';
            %             symbols = 2.^(0:log2(obj.M)-1)*data_2bit;
            %转换为十进制
            %             symbols=bi2de(data_2bit, 'left-msb');
            %             symbols(symbols == 0) = -3;
            %             symbols(symbols == 1) = -1;
            %             symbols(symbols == 2) = 1;
            %             symbols(symbols == 3) = 3;
        end

        function symbols_rand=rand_bits(obj)
            rng(obj.order);
            %参数：obj.prbsOrder，NSym，M
            symbols_rand=randi([0,1],log2(obj.M),obj.NSym);
            %             symbols_rand=randi([0 log2(obj.M)-1],obj.nModCarriers,obj.nPkts);
            %             symbols_rand = 2.^(0:log2(obj.M)-1)*data_bit;
        end


        function qam_signal=qam(obj,symbols)
            % 因为相位共轭关系，所以qam矩阵载波减半
            L=size(symbols,2)/2;
            qam_signal=qammod(symbols(:,1:L),obj.M,'InputType','bit','UnitAveragePower',1) ;
            qam_signal = reshape(qam_signal,obj.nModCarriers/2,[]);
        end

        %         function pcp_signal=phase_conj(obj,symbols)
        %             qam_signal=qammod(symbols,obj.M,'InputType','bit','UnitAveragePower',1) ;
        %             qam_signal = reshape(qam_signal,obj.nModCarriers,[]);
        %             % 所有信号的共轭
        %             pcp_signal=conj(qam_signal);
        %         end



        %ofdm信号生成
        function [ofdm_signal,postiveCarrierIndex] = ofdm(obj,qam_signal)
            % ofdm have two systerms as gapped OFDM and the interleaved OFDM
            % gapped OFDM :  1,2,3……size/2 装载数据
            % interleaved OFDM：1,3,5,7……size-1 装载数据
            %信号在矩阵中的位置
            postiveCarrierIndex = obj.nOffsetSub + (1:obj.nModCarriers);
            % nOffsetSub 行置零
            ofdmSig = ifft([zeros(obj.nOffsetSub,obj.nPkts);qam_signal;zeros(obj.fft_size-obj.nModCarriers-obj.nOffsetSub,obj.nPkts)]);
            % 第一行置零 ，第二行上载波信号，后续置零
            %             ofdmSig = ifft([zeros(1,obj.nPkts);qam_signal;zeros(obj.fft_size-obj.nModCarriers-1,obj.nPkts)]);
            ofdmSig = [ofdmSig(end-obj.nCP+1:end,:);ofdmSig];
            % 并串转换
            ofdmSig = ofdmSig(:);
            % 归一化
            scale_factor = max(max(abs(real(ofdmSig))),max(abs(imag(ofdmSig))));
            ofdm_signal = ofdmSig./scale_factor;
        end


        function     [ofdm_signal,postive_odd_index,postive_even_index,qam_mat]=ofdm_pcp_50(obj,qam_signal)

            % PC_pilot
            pcp_signal=conj(qam_signal);


            % 设置索引
            n_odd=1:2:obj.nModCarriers;
            n_even=2:2:obj.nModCarriers;
            % 得到新的QAM矩阵
            qam_mat=zeros(obj.nModCarriers,obj.nPkts);
            qam_mat(n_odd,:)=qam_signal;
            qam_mat(n_even,:)=pcp_signal;

            postive_odd_index=obj.nOffsetSub+n_odd;% 奇数放置Data
            postive_even_index=obj.nOffsetSub+n_even;%偶数放置phase data
            % nOffsetSub 行置零 ,OFDM_PCP信号
            X=zeros(obj.fft_size,obj.nPkts);
            X(postive_odd_index,:)=qam_signal;
            X(postive_even_index,:)=pcp_signal;
            % 转换为时域
            ofdmSig=ifft(X);
            ofdmSig = [ofdmSig(end-obj.nCP+1:end,:);ofdmSig];
            % 并串转换
            ofdmSig = ofdmSig(:);
            ofdm_signal=ofdmSig;
        end

        function     [ofdm_signal,index_data,index_pcp]=ofdm_pcp_percentage_33(obj,qam_signal)
            % 需要插入PCP的数量
            N=obj.nModCarriers*obj.percentage;
            % 只代表33% 的PCP
            n=2:1:N;
            index_pcp=obj.nOffsetSub+[2,3*n-1];
            % 找到需要进行conj的data
            index_data=index_pcp-1;


            % 初始信号
            X=[zeros(obj.nOffsetSub,obj.nPkts);qam_signal;zeros(obj.fft_size-obj.nModCarriers-obj.nOffsetSub,obj.nPkts)];
            % 进行提取,并进行转换
            X_pcp=conj(X(index_data));
            % 放置PCP
            X(index_pcp)=X_pcp;

            ofdmSig=ifft(X);
            % CP
            ofdmSig = [ofdmSig(end-obj.nCP+1:end,:);ofdmSig];
            % 并串转换
            ofdmSig = ofdmSig(:);

            % 归一化
            scale_factor = max(max(abs(real(ofdmSig))),max(abs(imag(ofdmSig))));
            ofdm_signal = ofdmSig./scale_factor;
        end

        function normalized_signal = nor(obj,signal)
            %归一化
            switch lower(obj.Type)
                %最大最小值归一化（Min-Max归一化）
                case 'max-min'
                    if isreal(signal)
                        min_val = min(signal);
                        max_val = max(signal);
                        normalized_signal = (signal - min_val) / (max_val - min_val);
                    else
                        real_min = min(real(signal));
                        real_max = max(real(signal));
                        imag_min = min(imag(signal));
                        imag_max = max(imag(signal));

                        % 实部和虚部的最大最小值归一化
                        normalized_real = (real(signal) - real_min) / (real_max - real_min);
                        normalized_imag = (imag(signal) - imag_min) / (imag_max - imag_min);

                        % 合并实部和虚部得到归一化后的复数信号
                        normalized_signal = normalized_real + 1i*normalized_imag;
                    end

                case 'max'
                    %除去最大值
                    if isreal(signal)
                        max_val = max(signal);
                        normalized_signal = signal/ max_val;
                    else
                        scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
                        normalized_signal = signal./scale_factor;
                    end

                case 'std'
                    %标准化
                    if isreal(signal)
                        s = std(signal);
                        m = mean(signal);
                        normalized_signal = (signal-m)/s;
                    else
                        real_mean = mean(real(signal));
                        real_std = std(real(signal));
                        imag_mean = mean(imag(signal));
                        imag_std = std(imag(signal));

                        % 实部和虚部的标准化
                        normalized_real = (real(signal) - real_mean) / real_std;
                        normalized_imag = (imag(signal) - imag_mean) / imag_std;

                        % 合并实部和虚部得到标准化后的复数信号
                        normalized_signal = normalized_real + 1i*normalized_imag;
                    end

                case 'unit'
                    %单位长度归一化
                    if isreal(signal)
                        % 计算实数信号的范数
                        norm_signal = norm(signal);

                        % 单位长度归一化
                        normalized_signal = signal / norm_signal;
                    else
                        % 计算每个复数的模长
                        abs_signal = abs(signal);

                        % 单位长度归一化
                        normalized_signal = signal ./ abs_signal;
                    end
            end
        end

        function hsqrt = get.hsqrt(obj)

            hsqrt = rcosdesign(obj.psfRollOff,obj.psfLength,obj.Nsam,obj.psfShape);  %设计根升余弦脉冲成型滤波器

        end

        function signal = MakeFilter(obj,qam_signal,hsqrt)
            % Upsample the data
            upam_signal = upsample(qam_signal,obj.Nsam);
            signal=conv(upam_signal,hsqrt,'same');
        end


        function psf = get.psf(obj)
            psf = fdesign.pulseshaping( obj.Nsam, obj.psfshape, ...
                'Nsym,beta', obj.psfLength, obj.psfRollOff);
        end

        function sig = DoFilter(obj,qam_signal)
            % Design the filter and get its frequencey response
            aux = design( obj.psf );
            hf = [ aux.Numerator, zeros( 1, ( obj.NSym - obj.psfLength ) * obj.Nsam - 1 ) ];
            hf = circshift( hf, [ 0, -obj.Nsam * obj.psfLength / 2 ] );
            % Upsample the data
            useq = upsample(qam_signal, obj.Nsam );
            sig = ifft( fft( useq ) .* fft( hf ) );
        end

        function plotSpectrum(obj)
            [y1,y2,Signal,~] = obj.Output();

            figure;
            f = obj.Fs * linspace( -0.5, 0.5, length(Signal) ).' / 1e9;
            plot( f, fftshift(10*log10(abs(fft(Signal))))); hold on;
            grid on;
            axis( [ obj.Fs / 2e9 * [ -1, 1 ], -70, 40 ] )

            figure;
            f = obj.Fs * linspace( -0.5, 0.5, length(y1) ).' / 1e9;
            plot( f, fftshift(10*log10(abs(fft(y1))))); hold on;
            grid on;
            axis( [ obj.Fs / 2e9 * [ -1, 1 ], -70, 40 ] )

            figure;
            f = obj.Fs * linspace( -0.5, 0.5, length(y2) ).' / 1e9;
            plot( f, fftshift(10*log10(abs(fft(y2))))); hold on;
            grid on;
            axis( [ obj.Fs / 2e9 * [ -1, 1 ], -70, 40 ] )
        end


    end
end