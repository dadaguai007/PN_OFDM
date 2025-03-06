% 复现:Dispersion-enhanced phase noise effects on reduced-guard-interval CO-OFDM transmission
clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')

% load data_kk_mat 加载接收矩阵（进行信道估计后）
% H向量大小为   符号数*1；

% OFDM 信号，行为carrier ，列为 symbols

M=16;
% first stage phase estimation
% 载波数
Num_Carrier=nn.nModCarriers;
% 选取合适导频进行phase 估计
pilotIndex= 1:2:Num_Carrier;
% 进行相位估计，并进行补偿
if strcmp(CPE_Status,'on')
    phase_compensation;
end

% 第一阶段后的接收信号
R=data_kk;
% 硬判决
% 提取每个载波的所有符号,进行硬判决

for index=1:size(R,1)
% 此数进行Weight_Decsion替换

    R_hat(index,:)=hard_decision(M,R(index,:));

end

% 共轭  conj(R_hat)

H=sum( R.* conj(R_hat));

% second stage phase estimation
phi= atan(imag(H)./real(H));

% 补偿
X=R.*...
    repmat(exp(-1j.*phi),size(R,1),1);



