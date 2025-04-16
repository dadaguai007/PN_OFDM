% 多载波脉冲成型

clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Project\Base_Code\Base\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);
SpS=nn.Fs/nn.Rs;
% 采样率
fs=nn.Fs;

% % 归一化
% scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
% signal = signal./scale_factor;


% 基本参数
Rs = 1.44e9;        % 符号速率 (1.44 Gbit/s)
Ts = 1/Rs;          % 符号周期 (0.694 ns)
alpha = 0.35;       % BTRC滚降因子
zeta = 10.47;       % Sinc-Lorentzian参数
gamma = zeta*Ts;    % Sinc-Lorentzian参数计算

% 时间向量
span = 10;          % 脉冲截断范围 (±10Ts)
t = (-span*Ts : Ts/100 : span*Ts)'; % 采样间隔为Ts/100

%% Rc

p_RC =  (1/Ts)*sinc(t/Ts) .* cos(pi*alpha*t/Ts) ./ (1-4*alpha^2*t.^2/Ts.^2);
% 归一化
% p_RC = p_RC /sqrt(sum(p_RC.^2));

% matlab的效果是一致的
filterCoeffs = rcFilterTaps(t, alpha, Ts);
%     filterCoeffs = filterCoeffs / sqrt(sum(filterCoeffs.^2));

%% BTRC脉冲生成
B = 1/(2*Ts);       % Nyquist带宽
beta = log(2)/(alpha*B); % BTRC参数

% BTRC时域表达式 (式4)
numerator = 4*beta*pi*t.*sin(2*pi*B*alpha*t) + 2*beta^2*cos(2*pi*B*alpha*t) - beta^2;
denominator = (4*pi^2*t.^2 + beta^2);
p_BTRC = 2*B * sinc(2*B*t) .* numerator ./ denominator;

% 归一化
p_BTRC = p_BTRC / max(abs(p_BTRC));

%% Sinc-Lorentzian脉冲生成 (式5-6)
alpha_minus = 1 - alpha; % 论文中设置为0

% 计算E(t) (式6)
term1 = alpha_minus * sinc(alpha_minus*t/Ts) / Ts;

term2_numerator = cos(pi*alpha_minus*t/Ts) - (pi*t/gamma).*sin(pi*alpha_minus*t/Ts);
term2_denominator = gamma + (pi^2*t.^2)/gamma;
term2 = term2_numerator ./ (term2_denominator + eps); % 避免除以0

E_t = term1 + term2;

% 合成脉冲 (式5)
p_SL = (gamma*Ts) / (gamma*alpha_minus + Ts) * sinc(t/Ts) .* E_t;

% 归一化
p_SL = p_SL / max(abs(p_SL));

%% 时域波形绘制
figure('Name','Time Domain');
hold on;
plot(t/Ts, p_BTRC, 'LineWidth',1.5);
title('BTRC Pulse (α=0.35)');
xlabel('Normalized Time (t/T_s)');
ylabel('Amplitude');
grid on;
xlim([-3 3]);

plot(t/Ts, p_SL, 'r', 'LineWidth',1.5);
% title('Sinc-Lorentzian Pulse (ζ=10.47)');
plot(t/Ts, p_RC, 'g', 'LineWidth',1.5);
grid on;
xlim([-3 3]);

%% 频谱分析
NFFT = 2^18;        % FFT点数
fs = 1/(t(2)-t(1)); % 采样频率

% BTRC频谱
P_BTRC = fftshift(fft(p_BTRC, NFFT));
f_BTRC = (-NFFT/2:NFFT/2-1)*(fs/NFFT)/1e9; % 频率单位GHz

% Sinc-Lorentzian频谱
P_SL = fftshift(fft(p_SL, NFFT));
P_SL_dB = 20*log10(abs(P_SL)/max(abs(P_SL)));

%Rc
P_RC = fftshift(fft(p_RC, NFFT));
P_RC_dB = 20*log10(abs(P_RC)/max(abs(P_RC)));

figure('Name','Spectrum');
plot(f_BTRC, 20*log10(abs(P_BTRC)/max(abs(P_BTRC))), 'LineWidth',1.5);
hold on;
plot(f_BTRC, P_SL_dB, 'r', 'LineWidth',1.5);
plot(f_BTRC, P_RC_dB, 'g', 'LineWidth',1.5);

title('Normalized Spectrum');
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
legend('BTRC','Sinc-Lorentzian','RC');
grid on;
xlim([-5 5]);