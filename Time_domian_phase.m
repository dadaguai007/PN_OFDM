% 完全时域迭代消除
% 复现：Flexible adjacent channel interference and phase noise suppression in energy-efficient OFDMA receivers
%% OFDM系统参数设置
clear; clc;close all;
N = 1024;               % 子载波总数
N_active = 600;         % 激活子载波数（中间300+两侧各150）
N_guard = (N - N_active)/2; % 保护子载波数
CP_len = 63;            % 循环前缀长度
mod_order = 16;         % 16QAM调制
beta = 350;             % 相位噪声3dB带宽（Hz）
Ts = 1/(15e3*N);        % 采样间隔（假设子载波间隔15kHz）
iter_num = 3;           % 迭代次数

%% 生成OFDM信号
data = randi([0 mod_order-1], N_active, 1);
tx_symbols = qammod(data, mod_order, 'UnitAveragePower', true); % 16QAM调制

% 子载波映射（中间激活，两侧保护）
tx_freq = [zeros(N_guard,1); tx_symbols; zeros(N_guard,1)];

% OFDM调制
tx_time = ifft(tx_freq, N);
tx_cp = [tx_time(end-CP_len+1:end); tx_time]; % 添加循环前缀

%% 信道模型（简化多径信道）
h = [1 0.5 0.3]; % 简单多径信道冲激响应
rx_channel = filter(h, 1, tx_cp);

%% 相位噪声生成（自由振荡器模型）
c = 4*pi*beta;          % 扩散率
phi = cumsum(sqrt(c*Ts)*randn(length(rx_channel),1)); % 维纳过程
rx_pn = rx_channel .* exp(1j*phi); % 添加相位噪声

%% 接收信号预处理
rx_signal = rx_pn(CP_len+1:end); % 去除循环前缀
rx_freq = fft(rx_signal, N);     % OFDM解调

%% 初始信道估计（假设完美已知）
H_est = fft(h, N).'; % 理想信道响应

%% 相位噪声消除算法核心
corrected_signal = rx_signal;
for iter = 1:iter_num
    % 信道均衡
    eq_freq = rx_freq ./ H_est;

    % 符号检测（硬判决）
    detected_symbols = qamdemod(eq_freq, mod_order, 'UnitAveragePower', true);
    tx_est = qammod(detected_symbols, mod_order, 'UnitAveragePower', true);

    % 重构时域信号
    tx_est_time = ifft([zeros(N_guard,1); tx_est(N_guard+1:end-N_guard); zeros(N_guard,1)], N);

    % 相位差估计
    phi_diff = angle(conj(tx_est_time) .* corrected_signal);

    % 低通滤波（移动平均滤波）
    LPF = ones(10,1)/10; % 示例：10阶FIR滤波器
    phi_est = filter(LPF, 1, phi_diff);

    % 相位补偿
    corrected_signal = corrected_signal .* exp(-1j*phi_est);

    % 更新频域信号
    rx_freq = fft(corrected_signal, N);
end

%% 邻道干扰处理（示例：频域滤波）
rx_freq_filtered = rx_freq;
rx_freq_filtered(1:N_guard) = 0;        % 滤除左侧邻道
rx_freq_filtered(end-N_guard+1:end) = 0;% 滤除右侧邻道

%% 最终符号检测
final_eq = rx_freq_filtered(N_guard+1:end-N_guard) ./ H_est(N_guard+1:end-N_guard);
received_data = qamdemod(final_eq, mod_order, 'UnitAveragePower', true);

%% 性能评估
SER = sum(data ~= received_data)/length(data);
disp(['Symbol Error Rate: ', num2str(SER)]);

%% 可视化
figure;
subplot(211); plot(real(tx_symbols));title('发射符号（实部）');
subplot(212); plot(real(final_eq), imag(final_eq), '.'); 
title('接收星座图'); axis square;