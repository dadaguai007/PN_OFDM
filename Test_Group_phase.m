% 复现:Dispersion-enhanced phase noise effects on reduced-guard-interval CO-OFDM transmission
clear;close all;clc;
addpath('Fncs\')
<<<<<<< HEAD
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
=======
addpath('D:\PhD\Codebase\')

% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
>>>>>>> 33b294655ad70c70f6633c1a87f6051c231dc07f
% load data_kk_mat 加载接收矩阵（进行信道估计后）
load OFDM_700km_fs32.mat
% H向量大小为   符号数*1；

% 原始信号
qam_signal_mat=qam_signal;
% H向量大小为   符号数*1；

% OFDM 信号，行为carrier ，列为 symbols

M=16;
% first stage phase estimation
% 载波数
Num_Carrier=size(data_kk_mat,1);
% 选取合适导频进行phase 估计
pilotIndex= 1:10:Num_Carrier;

% 数据采集
data_kk=data_kk_mat;

%%-----------------------------------  进行相位估计，并进行补偿   ----------------------------------------%%

% phase_compensation;
% phase_com_total_subcarrier;
% 归一化
data_rec=data_kk(:);
data_rec = data_rec./sqrt(mean(abs(data_rec(:)).^2));
% 第一阶段后的接收信号
R=reshape(data_rec,size(data_kk,1),[]);
% %%---------------------------------------  硬判决 ---------------------------------------------------------%%
% % 提取每个载波的所有符号,进行硬判决
% for index=1:size(R,1)
%     % 硬判
%     %     R_hat(index,:)=hard_decision_qam(M,R(index,:));
%     % Weight_Decision
%     R_hat(index,:)=Weighted_Decision(R(index,:));
% end

%%----------------------------------------- 分组  ------------------------------------------------------------%%
<<<<<<< HEAD
Group_Num = 10;
=======
Group_Num = 4;
>>>>>>> 33b294655ad70c70f6633c1a87f6051c231dc07f

for m=1:Num_Carrier/Group_Num
    % 每组数据的索引
    Num=(m-1)*Group_Num+1:1:m*Group_Num;
    phi=mean(R(Num,:)./qam_signal_mat(Num,:),1);
    phi_mean=angle(phi);

    X(Num,:)=R(Num,:).*...
        repmat(exp(-1j.*phi_mean),size(R(Num,:),1),1);

end




% 解码
OFDM_Decode_fun(X,ref_seq,M)