addpath(genpath('Fcns'));

nn=OFDMQAMN();
nn.DataType='rand';%两个选项：prbs，rand
nn.fft_size = 1024;
nn.nPkts = 1000;
nn.nCP = 32;
nn.nModCarriers = 250;
nn.nOffsetSub =250; 
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
 
f_idex=((nn.nModCarriers+nn.nOffsetSub)/nn.fft_size)*nn.Fs; 

% OFDM：1024*nPkts+32*nPkts