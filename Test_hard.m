clc;clear;
% 定义16QAM星座图
constellation = [-3-3i, -3-1i, -3+1i, -3+3i, ...
                 -1-3i, -1-1i, -1+1i, -1+3i, ...
                  1-3i,  1-1i,  1+1i,  1+3i, ...
                  3-3i,  3-1i,  3+1i,  3+3i];

% 生成一些示例的16QAM信号（这里简单假设是理想信号加上高斯噪声）
numSymbols = 1000; % 符号数量
noisePower = 0.1; % 噪声功率
txSymbols = randi([1, 16], 1, numSymbols); % 随机生成发送符号的索引
txSignal = constellation(txSymbols); % 生成发送信号
rxSignal = txSignal + sqrt(noisePower/2)*(randn(1, numSymbols) + 1i*randn(1, numSymbols)); % 接收信号，添加噪声

% 硬判决过程
rxSymbols = zeros(1, numSymbols); % 初始化判决后的符号索引
for i = 1:numSymbols
    distances = abs(rxSignal(i) - constellation).^2; % 计算接收信号点到所有星座点的距离
    [~, index] = min(distances); % 找到距离最近的星座点的索引
    rxSymbols(i) = index; % 记录判决结果
end

% 恢复出的信号
rxSignal_hardDecision = constellation(rxSymbols);

% 计算误符号率（SER）
errorCount = sum(rxSymbols ~= txSymbols);
SER = errorCount / numSymbols;
fprintf('误符号率 (SER): %.4f\n', SER);

% 绘制星座图
figure;
scatter(real(rxSignal), imag(rxSignal), 'b', 'filled', 'DisplayName', '接收信号');
hold on;
scatter(real(rxSignal_hardDecision), imag(rxSignal_hardDecision), 'r', 'filled', 'DisplayName', '判决后信号');
scatter(real(constellation), imag(constellation), 'ko', 'DisplayName', '星座点');
hold off;
xlabel('实部');
ylabel('虚部');
title('16QAM信号硬判决结果');
legend;
grid on;