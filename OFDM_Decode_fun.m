function OFDM_Decode_fun(data_kk_mat,ref_seq,M)

%归一化
data_kk=data_kk_mat(:);
data_kk = data_kk./sqrt(mean(abs(data_kk(:)).^2));

% 解码，计算误码
% 参考序列
ref_seq_1 =qamdemod(ref_seq,M,'OutputType','bit','UnitAveragePower',1);
ref_seq_1=ref_seq_1(:);
% 接收序列
yyy = data_kk;

yyy_1 = qamdemod(yyy,M,'OutputType','bit','UnitAveragePower',1);
yyy_1=yyy_1(:);
% 全部的序列进行解码
[ber1(squ_num),num1(squ_num),error_location] = CalcBER(yyy_1,ref_seq_1); %计算误码率
fprintf('Num of Errors = %d, BER = %1.7f\n',num1(squ_num),ber1(squ_num));
% 按符号解码，每个符号上的错码数
Calc_BER_mat;


end