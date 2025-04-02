%% LI-TD相位噪声抑制算法实现
% 论文: A time-domain phase noise mitigation algorithm for OFDM systems in the wireless backhaul links
% 作者: Peyman Neshaastegaran, Ming Jian
% 实现功能: 适用于高SNR无线回传链路的OFDM系统相位噪声抑制

% 代码作用： 如何初始化  估计相噪
clear all; close all; clc;

%% 系统参数设置
N = 4096;               % OFDM子载波数
N_A = 3300;             % 激活子载波数(根据5G-NR 400MHz带宽)
cp_len = 288;           % 循环前缀长度
M = 256;                % QAM调制阶数(256/1024)
SNR_dB = 30;            % 信噪比(dB)
PN_level = -85;         % 相位噪声水平(dBc/Hz @100kHz)
q = 20;                 % 滑动平均窗口大小
max_iter = 8;           % 最大迭代次数
pilot_ratio = 0.02;     % 导频比例(2%)

%% 生成相位噪声(维纳过程模型)
function phi = generate_phase_noise(N, PN_level)
    % 输入: N - 样本数, PN_level - 相位噪声水平(dBc/Hz)
    % 输出: phi - 相位噪声向量(弧度)
    
    beta = 10^(PN_level/10);    % 转换为线性值
    delta_phi = sqrt(2*pi*beta)*randn(1,N); % 相位增量
    phi = cumsum(delta_phi);    % 累积相位噪声
end

%% 初始TD相位噪声估计(基于DCT基)
function phi_est = initial_PN_estimation(y, S_pilot, pilot_pos, N, d)
    % 输入: y - 接收信号, S_pilot - 导频符号, pilot_pos - 导频位置
    %       N - FFT大小, d - DCT基向量数
    % 输出: phi_est - 初始相位噪声估计
    
    % 生成DCT基矩阵
    V = zeros(N,d);
    for i = 1:d
        V(:,i) = cos(pi*i/N*(0:N-1)');
    end
    
    % 计算矩阵M = F*diag(y)*V
    D_y = diag(y);
    M = fft(D_y * V, N);
    
    % 提取导频位置对应的行
    M_p = M(pilot_pos,:);
    S_p = S_pilot(:);
    
    % 计算系数估计(伪逆)
    gamma_est = pinv(M_p) * S_p;
    
    % 重构相位噪声
    phi_est = V * gamma_est;
end

%% LI-TD主算法
function [Y_est, phi_final] = LI_TD_PN_mitigation(y, H, N, M, q, max_iter, pilot_ratio)
    % 输入: y - 时域接收信号, H - 频域信道响应, N - FFT大小
    %       M - QAM阶数, q - 滑动窗口大小, max_iter - 最大迭代次数
    %       pilot_ratio - 导频比例
    % 输出: Y_est - 最终频域信号, phi_final - 最终相位噪声估计
    
    % 1. 初始参数设置
    N_pilot = round(N * pilot_ratio);    % 导频数量
    pilot_pos = randperm(N, N_pilot);    % 随机导频位置
    S_pilot = qammod(randi([0 M-1], N_pilot, 1), M, 'UnitAveragePower', true); % 导频符号
    
    % 2. 初始相位噪声估计(迭代0)
    d = 8;  % DCT基向量数
    phi_est = initial_PN_estimation(y, S_pilot, pilot_pos, N, d);
    
    % 3. 迭代处理
    for iter = 1:max_iter
        % 相位噪声补偿
        y_comp = y .* exp(-1j*phi_est);
        
        % OFDM解调
        Y_est = fft(y_comp, N);
        
        % 硬判决(使用最近邻解码)
        S_est = qamdemod(Y_est, M, 'UnitAveragePower', true);
        S_est_mod = qammod(S_est, M, 'UnitAveragePower', true);
        
        % 检查停止条件(连续两次判决相同)
        if iter > 1 && isequal(S_est, S_est_prev)
            break;
        end
        S_est_prev = S_est;
        
        % 信号重构
        r_est = ifft(H .* S_est_mod);
        
        % 计算z_k = y_k / r_est_k (公式7)
        z = y ./ r_est;
        
        % 滑动平均滤波
        num_blocks = ceil(N/q);
        phi_avg = zeros(num_blocks,1);
        for t = 1:num_blocks
            start_idx = (t-1)*q + 1;
            end_idx = min(t*q, N);
            z_block = z(start_idx:end_idx);
            phi_avg(t) = angle(mean(z_block)); % 公式11
        end
        
        % 线性插值
        x_blocks = linspace(1, N, num_blocks);
        x_samples = 1:N;
        phi_est = interp1(x_blocks, phi_avg, x_samples, 'linear');
    end
    
    phi_final = phi_est;
end

%% 主程序流程
% 1. 生成发送信号
data = randi([0 M-1], N_A, 1);
S = qammod(data, M, 'UnitAveragePower', true);

% 2. OFDM调制
X = ifft(S, N);

% 3. 添加循环前缀
X_cp = [X(end-cp_len+1:end); X];

% 4. 信道模型(简单AWGN信道)
H = ones(N,1); % 假设理想信道
h = ifft(H);

% 5. 生成相位噪声
phi_tx = generate_phase_noise(length(X_cp), PN_level);
phi_rx = generate_phase_noise(length(X_cp), PN_level);
phi_total = phi_tx + phi_rx;

% 6. 通过信道
y_channel = conv(X_cp, h);
y_channel = y_channel(1:length(X_cp)); % 保持长度不变

% 7. 添加相位噪声和AWGN
y_noisy = y_channel .* exp(1j*phi_total);
signal_power = mean(abs(y_noisy).^2);
noise_power = signal_power / (10^(SNR_dB/10));
noise = sqrt(noise_power/2) * (randn(size(y_noisy)) + 1j*randn(size(y_noisy)));
y_noisy = y_noisy + noise;

% 8. 去除循环前缀
y_noisy_cp = y_noisy(cp_len+1:end);

% 9. 执行LI-TD相位噪声抑制
[Y_est, phi_est] = LI_TD_PN_mitigation(y_noisy_cp, H, N, M, q, max_iter, pilot_ratio);

% 10. 性能评估
Y_ideal = fft(ifft(S).*H); % 理想无噪声无PN情况
MSE_PN = 10*log10(mean(abs(phi_est - phi_total(cp_len+1:end)').^2));
SNR_out = 10*log10(mean(abs(S).^2) / mean(abs(Y_est(1:N_A) - S).^2));

disp(['相位噪声估计MSE: ', num2str(MSE_PN), ' dB']);
disp(['输出SNR: ', num2str(SNR_out), ' dB']);

%% 可视化结果
figure;
subplot(2,1,1);
plot(phi_total(cp_len+1:end), 'b'); hold on;
plot(phi_est, 'r--');
legend('真实相位噪声', '估计相位噪声');
title('相位噪声估计对比');
xlabel('样本索引'); ylabel('相位(弧度)');

subplot(2,1,2);
plot(real(S), imag(S), 'bo'); hold on;
plot(real(Y_est(1:N_A)), imag(Y_est(1:N_A)), 'rx');
legend('发送符号', '接收符号');
title('星座图对比');
axis square;
