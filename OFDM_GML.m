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
    R_hat(index,:)=hard_decision(M,R(index,:));
end

% 分组
Group_Num = 30;

for m=1:Num_Carrier/Group_Num
    % 每组数据的索引
    Num=(m-1)*Group_Num+1:1:m*Group_Num;
    % 每组的H矩阵
    H(m,:)=sum( R(Num,:).* conj(R_hat(Num,:)));
    % second stage phase estimation
    phi(m,:)= atan(imag(H)./real(H));

end


% 补偿
for m=1:Num_Carrier/Group_Num
    % 每组数据的索引
    Num=(m-1)*Group_Num+1:1:m*Group_Num;
    % 每组之间进行补偿
    X(Num,:)=R(Num,:).*...
             repmat(exp(-1j.*phi(m,:)),Group_Num,1);
end
