% 复现:Dispersion-enhanced phase noise effects on reduced-guard-interval CO-OFDM transmission

clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')

    % 首先进行相位估计，并进行补偿
    if strcmp(CPE_Status,'on')
        phase_compensation;
    end

    R_hat=hard_decision(M,R);
H=R* R_hat';
phi= atan(imag(H)./real(H));


    % 进行判决
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

