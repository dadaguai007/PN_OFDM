


% pilotIndex=1:1:W;
phi_mean=angle((data_kk(pilotIndex,:)./...
    qam_signal_mat(pilotIndex,:)));


data_kk=data_kk.*...
    (exp(-1j.*phi_mean));