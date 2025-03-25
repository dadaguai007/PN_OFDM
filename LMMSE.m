% MATLAB 代码：基于线性最小均方误差（LMMSE）估计方法的相噪估计

% 清空工作区和命令窗口
clear;
clc;

% 系统参数设置
N = 64; % 子载波数量
Np = 8; % 导频数量
T = 1e-6; % 符号持续时间
Ts = T / N; % 采样间隔
R = 1e6; % 传输数据率
beta = 50; % 相噪线宽
SNR = 20; % 信噪比（dB）
numSymbols = 100; % 符号数量

% 相噪生成
sigma_u = sqrt(2 * pi * beta * T / N); % 相噪方差
phi_Tx = cumsum(randn(1, numSymbols * N) * sigma_u); % 发射机相噪
phi_Rx = cumsum(randn(1, numSymbols * N) * sigma_u); % 接收机相噪
phi = phi_Tx + phi_Rx; % 总相噪

% 信道模型（瑞利衰落信道）
h = sqrt(0.5) * (randn(N, numSymbols) + 1i * randn(N, numSymbols)); % 归一化信道系数

% 发射信号生成
X = randn(N, numSymbols) + 1i * randn(N, numSymbols); % 随机数据符号
X_pilot = randn(Np, numSymbols) + 1i * randn(Np, numSymbols); % 导频符号

% 接收信号生成
Y = zeros(N, numSymbols);
for m = 1:numSymbols
    % 相噪影响
    phi_m = phi((m-1)*N+1:m*N);
    X_m = ifft(X(:, m));
    X_m_noisy = X_m .* exp(1i * phi_m);
    X_m_noisy = fft(X_m_noisy);

    % 信道影响
    X_m_noisy = X_m_noisy .* h(:, m);

    % 加性高斯噪声
    noise = (randn(N, 1) + 1i * randn(N, 1)) / sqrt(2);
    noise = noise * (10^(-SNR/20));
    Y(:, m) = X_m_noisy + noise;
end

% LMMSE 估计
sigma_z = 10^(-SNR/20); % 噪声方差
R_p = exp(-pi * beta * abs(repmat(0:N-1, N, 1) - repmat((0:N-1)', N, 1)) * Ts); % 相噪相关矩阵

P_hat = zeros(N, numSymbols);
for m = 1:numSymbols
    % 构造接收信号矩阵
    Y_vec = Y(:, m);
    W = diag(h(:, m) .* X(:, m));

    % 计算 LMMSE 估计
    R_yy = W * R_p * W' + sigma_z^2 * eye(N);
    R_py = R_p * W';
    P_hat(:, m) = R_py / R_yy * Y_vec;
end

% 估计结果分析
phi_estimated = zeros(1, numSymbols);
for m = 1:numSymbols
    phi_estimated(m) = angle(P_hat(1, m));
end

% 绘制结果
figure;
subplot(2, 1, 1);
plot(1:numSymbols, angle(P_hat(1, :)), 'b-', 'LineWidth', 1.5);
title('估计的相噪');
xlabel('符号索引');
ylabel('相位（弧度）');
grid on;

subplot(2, 1, 2);
plot(1:numSymbols, phi(1:numSymbols), 'r--', 'LineWidth', 1.5);
hold on;
plot(1:numSymbols, phi_estimated, 'b-', 'LineWidth', 1.5);
legend('真实相噪', '估计相噪');
xlabel('符号索引');
ylabel('相位（弧度）');
grid on;