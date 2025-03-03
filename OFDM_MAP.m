% 复现:Low complexity blind detection in OFDM systems with phase noise

clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')

% 算法分为两个迭代部分：
% 相位估计迭代
% 接收信号进行估计迭代

R; % 接收信号

prior_mu; %相位噪声先验均值(向量)
prior_Sigma;%相位噪声先验协方差矩阵
sigma_sq; % 噪声方差

% % 相位噪声先验参数（假设零均值、单位协方差）
% prior_mu = zeros(4*M+2, 1);
% prior_Sigma = eye(4*M+2);



[X_prev,H] = initial_data_estimation(H, R) ;


% 构造T矩阵，由B矩阵为基本元素！！！
% 构造矩阵B (频域循环卷积的子矩阵)
B = construct_B(H, X_prev, Omega, N);

% 构造实虚合并矩阵T = [real(B); imag(B)]
T = [real(B); imag(B)];

% 转换为实值问题，接收信号分解为实部和虚部
R_real = [real(R); imag(R)];
% 计算相位噪声估计beta_hat（公式21） T'确实是共轭
term1 = real(T' * T) + (sigma_sq/2) * inv(prior_Sigma);
term2 = real(T' * R) + (sigma_sq/2) * inv(prior_Sigma) * prior_mu;

% 求向量bete_hat(估计得到)
beta_hat = term1 \ term2;

% 提取复数形式的相位噪声估计alpha_Omega（公式22）
alpha_Omega = beta_hat(1:2*M+1) + 1j * beta_hat(2*M+2:end);

% 步骤3: 数据符号更新
% 相位噪声补偿（公式27）
X_LS = compensate_phase_noise(R, H, alpha_Omega, Omega, N);

% 硬判决（公式28）
X_current = hard_decision(X_LS);




%% 子函数1: 初始数据估计（忽略相位噪声）
function [X_init,H] = initial_data_estimation(data_kk_ofdm,qam_signal,HK)
% 简化的LS均衡 + 硬判决
% channel estimation
rxTrainSymbol = data_kk_ofdm(:,1:nTrainSym);
qam_signal_mat=repmat(qam_signal,1,HK);
refTrainSymbol = qam_signal_mat(:,1:nTrainSym);
Hf = mean(rxTrainSymbol./refTrainSymbol,2);
% channel equalization(LS均衡)
data_kk = data_kk_ofdm.*repmat(1./Hf,1,nn.nPkts*HK);
% 信道响应矩阵：
H=repmat(Hf,1,nn.nPkts*HK);
% 硬判决
X_init=hard_decision(data_kk);
end


function rxSignal=hard_decision(M,R)
% 硬判决
const = qammod([0:M-1],M);
% 归一化
const=pnorm(const);
for i = 1:length(R)
    distances = abs(R(i) - const).^2; % 计算接收信号点到所有星座点的距离
    [~, index] = min(distances); % 找到距离最近的星座点的索引
    rxSymbols(i) = index; % 记录判决结果
end
% 恢复出的信号
rxSignal= const(rxSymbols);
end
%% 子函数2: 构造矩阵B（频域循环卷积子矩阵）
function B = construct_B(R, Omega, N)
% 输入:
% H: 信道频域响应 (N x 1)
% X: 当前数据符号估计 (N x 1)
% Omega: 相位噪声频谱分量索引 (例如[N-1, 0, 1]对应循环索引)
% N: 子载波数


% 直接用接收的信号进行代替？
HX=R;
%     HX = H .* X;  % 逐元素相乘
num_terms = length(Omega);
B = zeros(N, num_terms);

for idx = 1:num_terms
    shift = Omega(idx);
    % 循环移位HX，对应频域循环卷积
    B(:, idx) = circshift(HX, shift);
end
end

%% 子函数3: 相位噪声补偿（频域操作）
function X_LS = compensate_phase_noise(R, H, alpha_Omega, Omega, N)
% 输入:
% R: 接收信号 (N x 1)
% H: 信道频域响应 (N x 1)
% alpha_Omega: 估计的相位噪声频谱分量 (2M+1 x 1)
% Omega: 相位噪声索引集合
% N: 子载波数

% 确认输出信号大小！！！

X_LS = zeros(N, 1);
%     M = (length(Omega) - 1)/2;

for k = 0:N-1
    sum_term = 0;

    for l_idx = 1:length(Omega)
        l = Omega(l_idx);
        % 计算 (k-l) mod N
        idx = mod(k-l, N);
        % 取alpha的共轭---对应公式27中的alpha^*[(-l)]
        alpha_conj = conj(alpha_Omega(l_idx));
        sum_term = sum_term + alpha_conj * R(idx + 1);  % MATLAB索引从1开始
    end
    % 从接收的信号中去除相噪
    X_LS(k + 1) = sum_term / H(k+1);
end
end
