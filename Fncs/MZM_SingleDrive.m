function [Eout] = MZM_SingleDrive(Ein,VRF,Vbias,Vpi,gamma)
%UNTITLED8 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
VRF1 = VRF/2;
VRF2 = -VRF/2;
Vbias1 = Vbias/2;
Vbias2 = -Vbias/2;
Eout = MZM_DualDrive(Ein,VRF1,VRF2,Vbias1,Vbias2,Vpi,gamma);
end

