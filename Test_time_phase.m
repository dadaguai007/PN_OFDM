%% 系统参数设置
clear all; clc;

% 基本参数
N = 4096;           % OFDM总子载波数
L = 4;              % 子载波分组数
M = N/L;            % 每组子载波数
L_cp = 4;           % 频域循环前缀长度
L_cs = 4;           % 频域循环后缀长度
SCS = 120e3;        % 子载波间隔(Hz)
fc = 100e9;         % 载波频率(Hz)
modOrder = 4;       % 调制阶数(QPSK)
maxIter = 10;       % 最大迭代次数
errThresh = 1e-3;   % 迭代停止阈值

%% 生成发送信号
dataBits = randi([0 1], N*log2(modOrder), 1); % 生成随机比特流
modSymbols = qammod(dataBits, modOrder, 'InputType', 'bit', 'UnitAveragePower', true); % QPSK调制

%% 子载波分组及循环前缀/后缀添加
txSymbols = reshape(modSymbols, M, L); % 将符号分配到L个子载波组

% 添加频域循环前缀和后缀
extendedSymbols = zeros(M+L_cp+L_cs, L);
for l = 1:L
    group = txSymbols(:, l);
    % 添加循环前缀（复制末尾L_cp个符号）
    cp = group(end-L_cp+1:end); 
    % 添加循环后缀（复制开头L_cs个符号）
    cs = group(1:L_cs);       
    extendedSymbols(:, l) = [cp; group; cs];
end

%% OFDM调制
timeSignal = ifft(extendedSymbols, N); % IFFT变换
txSignal = timeSignal(:); % 串行化

%% 信道效应（添加相位噪声）
% 相位噪声生成（维纳过程模型）
theta = cumsum(0.1*randn(size(txSignal))); % 相位噪声随机游走
pn = exp(1j*theta);        % 相位噪声时域表达式

% 接收信号（考虑相位噪声影响）
rxSignal = txSignal .* pn; % 乘性相位噪声
rxSignal = awgn(rxSignal, 20); % 添加高斯白噪声

%% 接收端处理
% OFDM解调
rxSymbols = reshape(rxSignal, N, L);
rxFreq = fft(rxSymbols); % FFT变换

% 移除循环前缀/后缀
processedGroups = zeros(M, L);
for l = 1:L
    fullGroup = rxFreq(:, l);
    processedGroups(:, l) = fullGroup(L_cp+1:end-L_cs); % 去除CP/CS
end

%% 相位噪声抑制(SIC算法)
decodedSymbols = zeros(M, L);
for l = 1:L
    currentGroup = processedGroups(:, l);
    prevEst = currentGroup;  % 初始化前次估计
    iter = 1;                % 迭代计数器
    error = inf;             % 初始化误差
    
    while iter <= maxIter && error > errThresh
        % 硬判决
        hardDecision = qamdemod(prevEst, modOrder, 'OutputType', 'approx', 'UnitAveragePower', true);
        
        % 重构时域信号
        timeEst = ifft(hardDecision); % IFFT变换
        
        % 估计相位噪声
        theta_est = angle(timeEst ./ ifft(prevEst)); % 相位差估计
        pn_est = exp(1j*theta_est); % 相位噪声估计
        
        % 频域响应计算并截断
        pn_freq = fft(pn_est);
        pn_freq(L_cp+1:end-L_cs) = 0; % 保留主瓣分量
        
        % 干扰消除
        updatedEst = fft(ifft(currentGroup) .* conj(pn_est)); % 频域解卷积
        updatedEst = updatedEst(1:M); % 截取有效子载波
        
        % 计算误差
        error = mean(abs(updatedEst - prevEst).^2);
        prevEst = updatedEst; % 更新估计值
        iter = iter + 1;
    end
    decodedSymbols(:, l) = prevEst; % 存储最终估计
end

%% 信号解调
finalSymbols = decodedSymbols(:); % 合并子载波组
rxBits = qamdemod(finalSymbols, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);

%% 性能评估
ber = sum(rxBits ~= dataBits)/numel(dataBits);
disp(['误码率: ', num2str(ber)]);
