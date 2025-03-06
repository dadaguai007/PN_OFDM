% 复现:Dispersion-enhanced phase noise effects on reduced-guard-interval CO-OFDM transmission
clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
% load data_kk_mat 加载接收矩阵（进行信道估计后）
load OFDM_4800km.mat
load OFDM_4800km_label.mat
% H向量大小为   符号数*1；

% 原始信号
qam_signal_mat=qam_signal;

% OFDM 信号，行为carrier ，列为 symbols

M=16;
% first stage phase estimation
% 载波数
Num_Carrier=size(data_kk_mat,1);
% 选取合适导频进行phase 估计
pilotIndex= 1:2:Num_Carrier;

% 数据采集
data_kk=data_kk_mat;
% 进行相位估计，并进行补偿

phase_compensation;


% 第一阶段后的接收信号
R=data_kk;
% 硬判决
% 提取每个载波的所有符号,进行硬判决

for index=1:size(R,1)

    % 硬判决
%     R_hat(index,:)=hard_decision_qam(M,R(index,:));
    % Weight_Decision
    R_hat(index,:)=Weighted_Decision(R(index,:));
end

% 共轭  conj(R_hat)

H=sum( R.* conj(R_hat));

% second stage phase estimation
phi= atan(imag(H)./real(H));

% 补偿
X=R.*...
    repmat(exp(-1j.*phi),size(R,1),1);


% 解码
OFDM_Decode_fun(X,ref_seq,M)
