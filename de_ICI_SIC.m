% 载波分组，时域相噪迭代消除
% 复现：Mitigation of phase noise-induced ICI at THz bands using CP-OFDM PT-RS signals
%% 参数设置
clear; close all; clc;

% 系统参数
N = 4096;           % 总子载波数
M = 256;            % 每个子载波组的长度
L_CP = 4;           % FD-CP长度
L_CS = 4;           % FD-CS长度
L = N / M;          % 子载波组数
SCS = 120e3;        % 子载波间距 (Hz)
fc = 100e9;         % 载波频率 (100 GHz)
modOrder = 16;      % 调制阶数 (16QAM)
EbN0_dB = 20;       % Eb/N0 (dB)
numSymbols = 10;    % OFDM符号数
maxIter = 10;       % 最大迭代次数
epsilon_T = 1e-3;   % 迭代终止阈值

% 相位噪声参数
FoM = -110;         % 振荡器品质因数 (dBc/Hz)
f_offset = 1e6;     % 相位噪声偏移频率 (Hz)
PN_model = 'Model1'; % 相位噪声模型选择

%% 发射端处理
% 生成随机数据符号
dataBits = randi([0 1], N*log2(modOrder), numSymbols);
modSymbols = qammod(dataBits, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);

% 子载波分区与FD-CP/CS插入
txSymbols = zeros(N + (L_CP + L_CS)*L, numSymbols);
for symIdx = 1:numSymbols
    % 将数据划分为L个子载波组
    groupedSymbols = reshape(modSymbols(:, symIdx), M, L);

    % 对每个子载波组添加FD-CP和FD-CS
    for l = 1:L
        group = groupedSymbols(:, l);
        cp = group(end-L_CP+1:end);      % 取后L_CP个符号作为CP
        cs = group(1:L_CS);              % 取前L_CS个符号作为CS
        txGroup = [cp; group; cs];       % 组合CP+数据+CS
        txSymbols((l-1)*(M+L_CP+L_CS)+1 : l*(M+L_CP+L_CS), symIdx) = txGroup;
    end
end

% OFDM调制（IFFT）
txSignal = ifft(txSymbols, N, 1);

%% 信道模型（含相位噪声）
% 生成相位噪声
theta = generatePhaseNoise(N, numSymbols, FoM, fc, SCS, PN_model);

% 应用相位噪声和AWGN
rxSignal = zeros(size(txSignal));
for symIdx = 1:numSymbols
    % 相位噪声乘性效应
    corruptedSignal = txSignal(:, symIdx) .* exp(1j * theta(:, symIdx));

    % 添加高斯噪声
    SNR = EbN0_dB + 10*log10(log2(modOrder)) - 10*log10(N/(N + (L_CP+L_CS)*L));
    noisePower = 10^(-SNR/10);
    noise = sqrt(noisePower/2) * (randn(size(corruptedSignal)) + 1j*randn(size(corruptedSignal)));
    rxSignal(:, symIdx) = corruptedSignal + noise;
end

%% 接收端处理
% 去时域CP并FFT
rxSymbols = fft(rxSignal, N, 1);

% 信道均衡（理想假设）
eqSymbols = rxSymbols; % 实际中需使用LMMSE均衡

% CPE补偿（利用PT-RS，此处简化为理想补偿）
cpe = angle(mean(eqSymbols(1:M:end, :))); % 假设PT-RS在固定位置
eqSymbols = eqSymbols .* exp(-1j * cpe);

% 子载波重组与去FD-CP/CS
rxGroups = zeros(M, L, numSymbols);
for symIdx = 1:numSymbols
    fullSymbol = eqSymbols(:, symIdx);
    for l = 1:L
        startIdx = (l-1)*(M+L_CP+L_CS) + L_CP + 1;
        endIdx = startIdx + M - 1;
        rxGroups(:, l, symIdx) = fullSymbol(startIdx:endIdx);
    end
end

%% 迭代SIC去ICI算法
decodedSymbols = zeros(size(modSymbols));
for symIdx = 1:numSymbols
    for l = 1:L
        Xk = rxGroups(:, l, symIdx);
        P_est = zeros(M, 1); % 初始化相位噪声系数

        for iter = 1:maxIter
            % 硬判决
            X_hat = qamdemod(Xk, modOrder, 'OutputType', 'approx', 'UnitAveragePower', true);

            % 重建时域信号
            x_hat = ifft(X_hat);

            % 估计相位噪声
            theta_est = angle(rxGroups(:, l, symIdx) ./ x_hat);
            P_est = fft(theta_est);

            % 截断高频分量
            P_est(L_CP+1:end-L_CS) = 0;

            % 频域解卷积
            Xk_new = ifft(fft(Xk) ./ P_est);

            % 计算误差
            epsilon = norm(Xk_new - Xk) / sqrt(M);
            if epsilon < epsilon_T
                break;
            end
            Xk = Xk_new;
        end

        decodedSymbols((l-1)*M+1 : l*M, symIdx) = Xk;
    end
end

%% BER计算
% 解调
decodedBits = qamdemod(decodedSymbols(:), modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
ber = sum(decodedBits ~= dataBits(:)) / numel(dataBits);
disp(['BER: ', num2str(ber)]);

%% 相位噪声生成函数
function theta = generatePhaseNoise(N, numSymbols, FoM, fc, SCS, model)
    % 根据模型生成相位噪声
    Ts = 1/(N*SCS); % 符号持续时间
    t = (0:N-1)*Ts;
    theta = zeros(N, numSymbols);

    switch model
        case 'Model1' % 参考3GPP模型
            phaseNoisePSD = @(f) 10^(FoM/10) ./ (f.^2); % 相位噪声功率谱密度
            for symIdx = 1:numSymbols
                f = linspace(-SCS/2, SCS/2, N);
                S_phi = phaseNoisePSD(abs(f));
                phi = sqrt(S_phi) .* (randn(1,N) + 1j*randn(1,N));
                theta(:,symIdx) = real(ifft(phi));
            end
        otherwise
            error('Unsupported PN model');
    end
end


%% 修改代码

%% 迭代SIC去ICI算法（修改部分）
% decodedSymbols = zeros(size(modSymbols));
% for symIdx = 1:numSymbols
%     for l = 1:L
%         Xk = rxGroups(:, l, symIdx);
%         P_est = zeros(M, 1); % 初始化相位噪声系数
% 
%         for iter = 1:maxIter
%             % 硬判决
%             X_hat = qamdemod(Xk, modOrder, 'OutputType', 'approx', 'UnitAveragePower', true);
% 
%             % 重建时域信号
%             x_hat = ifft(X_hat);
% 
%             % 估计相位噪声
%             theta_est = angle(rxGroups(:, l, symIdx) ./ x_hat);
%             P_est = fft(theta_est);
% 
%             % --- 修改1：截断高频分量 ---
%             P_est_trunc = P_est;
%             P_est_trunc(L_CP+1 : end-L_CS) = 0;
% 
%             % --- 修改2：循环翻转（保留P0，翻转其余部分）---
%             P_est_flipped = [P_est_trunc(1); flip(P_est_trunc(2:end))];
% 
%             % --- 修改3：频域解卷积（使用翻转后的相位噪声响应）---
%             Xk_freq = fft(Xk);
%             Xk_new_freq = Xk_freq ./ (P_est_flipped + eps); % 添加正则化
% 
%             % 转换回时域
%             Xk_new = ifft(Xk_new_freq);
% 
%             % 计算误差
%             epsilon = norm(Xk_new - Xk) / sqrt(M);
%             if epsilon < epsilon_T
%                 break;
%             end
%             Xk = Xk_new;
%         end
% 
%         decodedSymbols((l-1)*M+1 : l*M, symIdx) = Xk;
%     end
% end