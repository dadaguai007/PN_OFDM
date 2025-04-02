% Dual pilot express the TX LPN

clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
% load data_kk_mat 加载接收矩阵（进行信道估计后）
load OFDM_700km_fs32.mat
% H向量大小为   符号数*1；

% 原始信号
qam_signal_mat=qam_signal;
% H向量大小为   符号数*1；

% OFDM 信号，行为carrier ，列为 symbols

M=16;
% 载波数
Num_Carrier=size(data_kk_mat,1);
% 选取合适导频进行phase 估计
pilotIndex= [1,10];

% 数据采集
data_kk=data_kk_mat;
%%-----------------------------------  进行相位估计，并进行补偿   ----------------------------------------%%
% 归一化
data_rec=data_kk(:);
data_rec = data_rec./sqrt(mean(abs(data_rec(:)).^2));
% 第一阶段后的接收信号
R=reshape(data_rec,size(data_kk,1),[]);


% Dual pilot phase
phi_k1=angle((data_kk(pilotIndex(1),:)./...
    qam_signal_mat(pilotIndex(1),:)));

phi_k2=angle((data_kk(pilotIndex(2),:)./...
    qam_signal_mat(pilotIndex(2),:)));
% 相位差
detal_phase1 = phi_k2-phi_k1;

c = 299792458;
c_kms = c/ 1e3 ;   % speed of light (vacuum) in km/s
Fc = 193.1e12;
lamba = c_kms / Fc ;
D=16;
beta2 = -(D * lamba.^2) / (2 * pi * c_kms);
T=1/32e9;
L=700;
detal_f=32e9/1023;
N=diff(pilotIndex)-1;
detal_w=2*pi*N*detal_f;
% 系数
E=beta2*L*detal_w/T;
% Tx LPN

for i=1:size(data_kk,2)

    index=1:i;
    LPN(i)=1/E*sum(detal_phase1(index));

end

% LPN 过大，方法应该错误
figure;
plot(LPN)

%  补延时，使用iqdelay函数