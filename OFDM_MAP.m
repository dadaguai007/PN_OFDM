% 复现:Low complexity blind detection in OFDM systems with phase noise

clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')

% 算法分为两个迭代部分：
% 相位估计迭代
% 接收信号进行估计迭代

% 接收的信号R 矩阵形状为 载波*符号

% 接下来的矩阵操作都是 符号* 载波

M=16;
% 假设接收矩阵为N*N
% ICI 收到2M+1 载波影响
R; 


% sigma_sq; % 噪声方差

% % % 相位噪声先验参数 大小设置完成
% prior_mu = zeros(4*M+2, 1); % 假设零均值 ， 
% % prior_Sigma = eye(4*M+2);


% 问题: 先验参数的估计：
% 参数设定
m = 1;                      % 相位噪声低频分量数
num_terms = 2*m + 1;        % Omega中的分量数
sigma_l = [0.1, 0.2, 0.05]; % 各频谱分量的标准差（示例）

% 先验均值（假设为零均值）%相位噪声先验均值(向量)
prior_mu = zeros(2*num_terms, 1);

% 协方差矩阵（假设频谱分量独立且Γ=0）
Sigma_alpha = diag(sigma_l.^2);      % α_Ω的协方差矩阵
Sigma_R = real(Sigma_alpha) / 2;     % 实部协方差
Sigma_I = real(Sigma_alpha) / 2;     % 虚部协方差
Sigma_RI = -imag(Sigma_alpha) / 2;   % 实虚交叉协方差
Sigma_IR = imag(Sigma_alpha) / 2;    % 虚实交叉协方差

% 构造Σβ %相位噪声先验协方差矩阵
prior_Sigma = [Sigma_R, Sigma_RI; 
              Sigma_IR, Sigma_I];

[X_prev,H] = initial_data_estimation(H, R) ;

for iter = 1:I_max
% 构造T矩阵，由B矩阵为基本元素！！！
% 构造矩阵B (频域循环卷积的子矩阵)
% B矩阵有两种方式： 全矩阵，以及2M+1载波的矩阵

B = construct_B(H, X_prev, Omega, N);

% 构造实虚合并矩阵 N*（4M+2）
T = [B; 1j*B];

% beta 应为 4M+2的长度

% 计算相位噪声估计beta_hat（公式21） T'确实是共轭
term1 = real(T' * T) + (sigma_sq/2) * inv(prior_Sigma);
term2 = real(T' * R) + (sigma_sq/2) * inv(prior_Sigma) * prior_mu;

% 求向量bete_hat(估计得到)
beta_hat = term1 \ term2;  % 等价于 inv(term1) * term2

% 提取复数形式的相位噪声估计alpha_Omega（公式22）
alpha_Omega = beta_hat(1:2*M+1) + 1j * beta_hat(2*M+2:end);

% 步骤3: 数据符号更新
% 相位噪声补偿（公式27）
X_LS = compensate_phase_noise(R, H, alpha_Omega, Omega, N);

% 硬判决（公式28）
X_current = hard_decision(X_LS);

% 检查终止条件：数据符号变化小于阈值（使用范数进行判断）
if norm(X_current - X_prev) < epsilon
    break;
end
X_prev = X_current;
end
X_hat = X_current;





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

% 此时的data_kk 是一个N*N的矩阵

% 信道响应矩阵： N* 1 的矩阵（按照符号长度定义）
H=repmat(Hf,1,nn.nPkts*HK);
% 硬判决
X_init=hard_decision(data_kk);
end



%% 子函数2: 构造矩阵B（频域循环卷积子矩阵）
% B矩阵应为N*2M+1

% 如何处理循环卷积矩阵
function B_index = construct_B(H,X, Omega, N,k)
% 输入:
% H: 信道频域响应 (N x 1)
% X: 当前数据符号估计 (N x 1)
% Omega: 相位噪声频谱分量索引 (例如[N-1, 0, 1]对应循环索引)
% N: 子载波数
% k是载波索引

% 取所有载波的第一个符号  →→→  后续可扩展为所有载波的单独一个符号
H=H(1,:);
X=X(1,:);
num_terms = length(Omega);
B = zeros(N, num_terms);
% ICI影响区间
M=(num_terms-1)/2;
%  影响的载波索引
Carrier_index=k-M:1:k+M;
% 只取2M+1的载波
% HX=R(:,Carrier_index);

% Hx为向量形式(应先形成矩阵，再从输出向量B中进行取值)
% HX = H(:,Carrier_index) .* X(:,Carrier_index);  % 逐元素相乘（刚估计出的X*信道矩阵）
% 生成信道的对角矩阵 以及接收信号的列向量 进行相乘，得到接收信号状态
HX=diag(H)*X.';

for idx = 1:num_terms
%     shift = Omega(idx);
    % 循环移位HX，对应频域循环卷积
    B(:, idx) = circshift(HX, idx-1);
end
% 取收到影响载波对应的B矩阵
B_index=B(:,Carrier_index);
end

%% 子函数3: 相位噪声补偿（频域操作）
function X_LS = compensate_phase_noise(R, H, alpha_Omega, Omega, N)
% 输入:
% R: 接收信号 (N x 1)
% H: 信道频域响应 (N x 1)
% alpha_Omega: 估计的相位噪声频谱分量 (2M+1 x 1)
% Omega: 相位噪声索引集合
% N: 子载波数


for k = 0:N-1
    sum_term = 0;

    for l_idx = 1:length(Omega)

        l = Omega(l_idx);

        % ICI影响区间
        M=(l-1)/2;
        %  影响的载波索引
        Carrier_index=k-M:1:k+M;

        % 取alpha的共轭---对应公式27中的alpha^*[(-l)]，循环矩阵的缘故
        alpha_conj = conj(alpha_Omega(:,l_idx));
        sum_term = sum_term + alpha_conj * R(:,Carrier_index(l_idex));  % 受到影响的载波索引
    end
    % 从接收的信号中去除相噪(对每一个载波单独计算)
    X_LS(:,k + 1) = sum_term ./ H(:,k+1);
end
end
