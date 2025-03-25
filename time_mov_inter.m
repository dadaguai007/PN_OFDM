clc; clear; close all;

%% 参数设置
N = 4096;               % FFT点数
N_active = 3300;        % 有效子载波数
cp_len = 288;           % 循环前缀长度
M = 256;                % QAM调制阶数（256或1024）
SNR_dB = 30;            % 信噪比（dB）
phase_noise_level = -85;% 相位噪声水平（dBc/Hz @100kHz）
q = 20;                 % 移动平均窗口大小
i_max = 8;              % 最大迭代次数

%% 生成OFDM信号
data = randi([0 M-1], N_active, 1);        % 生成随机数据
S = qammod(data, M, 'UnitAveragePower', true); % QAM调制
X = zeros(N, 1);
X(1:N_active/2) = S(1:N_active/2);         % 有效子载波分配
X(end-N_active/2+1:end) = S(N_active/2+1:end);
x = ifft(X, N);                            % IFFT
x_cp = [x(end-cp_len+1:end); x];           % 添加循环前缀

%% 模拟信道（添加相位噪声和AWGN）
% 生成维纳相位噪声
beta = 10^(phase_noise_level/10);          % 相位噪声方差
phi_tx = cumsum(sqrt(beta)*randn(size(x_cp))); % 发射端相位噪声
phi_rx = cumsum(sqrt(beta)*randn(size(x_cp))); % 接收端相位噪声
phi_total = phi_tx + phi_rx;               % 总相位噪声

% 添加相位噪声和AWGN
y_cp = x_cp .* exp(1j*phi_total);          % 时域信号受相位噪声影响
y_cp = awgn(y_cp, SNR_dB, 'measured');     % 添加AWGN

% 去除循环前缀
y = y_cp(cp_len+1:end);

%% 初始相位噪声估计（迭代0）
[phi_hat_initial, S_hat_initial] = initial_PN_estimation(y, N, N_active, M);

%% LI-TD迭代优化
phi_hat = phi_hat_initial;
S_hat = S_hat_initial;
for iter = 1:i_max
    [phi_hat, S_hat] = LI_TD_iteration(y, phi_hat, S_hat, N, N_active, M, q);

    % 检查停止条件（符号估计不再变化）
    if iter > 1 && isequal(S_hat, S_hat_prev)
        break;
    end
    S_hat_prev = S_hat;
end

%% 结果评估
Y_hat = fft(y .* exp(-1j*phi_hat), N);     % 相位补偿后的频域信号
Y_hat_active = [Y_hat(1:N_active/2); Y_hat(end-N_active/2+1:end)];
data_hat = qamdemod(Y_hat_active, M, 'UnitAveragePower', true);
SER = sum(data ~= data_hat) / N_active;     % 符号错误率
fprintf('符号错误率 (SER): %.4f\n', SER);

%% 绘制相位噪声估计误差
figure;
plot(phi_total(cp_len+1:end) - phi_hat, 'r');
title('相位噪声估计误差');
xlabel('样本索引'); ylabel('误差 (rad)');



function [phi_hat, S_hat] = initial_PN_estimation(y, N, N_active, M)
    % 使用DCT基向量进行初始估计（参考论文Section II-A）
    d = 8;                          % 基向量数量（与论文参数一致）
    V = dctmtx(N)';                 % 生成DCT基矩阵（N x d）
    V = V(:,1:d);                   % 选择前d个基向量

    % 构造观测矩阵（参考公式(5)）
    Y = diag(y);
    F = fft(eye(N))/sqrt(N);        % FFT矩阵
    Lambda = eye(N);                % 假设理想信道（可根据实际信道修改）
    M_matrix = F * Y * V;
    gamma_hat = pinv(M_matrix(1:N_active,:)) * S_hat_initial; % 伪逆运算

    % 估计相位噪声
    phi_hat = angle(V * gamma_hat); % 相位估计

    % 初始符号判决
    Y_hat_initial = fft(y .* exp(-1j*phi_hat), N);
    Y_hat_active = [Y_hat_initial(1:N_active/2); Y_hat_initial(end-N_active/2+1:end)];
    S_hat = qamdemod(Y_hat_active, M, 'UnitAveragePower', true);
end



function [phi_hat_new, S_hat_new] = LI_TD_iteration(y, phi_hat_prev, S_hat_prev, N, N_active, M, q)
    % 重构时域信号（公式(1)）
    X_hat_prev = zeros(N, 1);
    X_hat_prev(1:N_active/2) = S_hat_prev(1:N_active/2);
    X_hat_prev(end-N_active/2+1:end) = S_hat_prev(N_active/2+1:end);
    x_hat_prev = ifft(X_hat_prev, N);

    % 计算z_k（公式(7)-(10)）
    z = y ./ x_hat_prev;            % 近似为单载波模型

    % 移动平均滤波
    z_filtered = movmean(z, q);     % 移动平均窗口大小为q

    % 线性插值
    t = 1:q:N;                      % 分段点
    phi_segments = angle(z_filtered(t)); % 分段相位估计
    phi_hat_new = interp1(t, phi_segments, 1:N, 'linear'); % 线性插值

    % 更新符号估计
    Y_hat_new = fft(y .* exp(-1j*phi_hat_new), N);
    Y_hat_active = [Y_hat_new(1:N_active/2); Y_hat_new(end-N_active/2+1:end)];
    S_hat_new = qamdemod(Y_hat_active, M, 'UnitAveragePower', true);
end