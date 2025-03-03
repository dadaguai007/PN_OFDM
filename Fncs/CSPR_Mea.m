% CSPR值计算

% theory
P_carrier=A.^2;
P_signal=signalpower(Amp*signal_ofdm);
CSPR = P_carrier/(P_signal);
CSPR=10*log10(CSPR);
fprintf('CSPR1 = %1.7f\n',CSPR);

