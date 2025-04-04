addpath(genpath('Fcns'));
% 带有相位共轭项的OFDM信号
nn=OFDMQAMN_phase_conjugated();
nn.DataType='rand';%两个选项：prbs，rand
nn.fft_size = 1024;
nn.nPkts = 1000;
nn.nCP = 32;
nn.nModCarriers = 300;
nn.nOffsetSub =2; 
% total symbol
nn.NSym = nn.nModCarriers*nn.nPkts;
nn.order = 2;
nn.M = 16;
nn.prbsOrder = 15;
nn.Rs = 64e9;
nn.Fs = 64e9;
nn.Nsam =nn.Fs/nn.Rs ;
nn.psfRollOff=0.01;
nn.psfLength=256;
nn.psfShape='sqrt';
nn.psfshape='Raised Cosine';
 nn.len= (nn.fft_size+nn.nCP)*nn.nPkts; % OFDM 符号长度
nn.dataCarrierIndex=nn.nOffsetSub+(1:nn.nModCarriers);

f_idex=((nn.nModCarriers+nn.nOffsetSub)/nn.fft_size)*nn.Fs; 
% 载波位置
postiveCarrierIndex=nn.nOffsetSub+1:1:nn.nModCarriers+nn.nOffsetSub;

% OFDM：1024*nPkts+32*nPkts