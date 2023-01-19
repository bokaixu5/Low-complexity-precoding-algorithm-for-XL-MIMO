close all;
clearvars;
j=sqrt(-1);
S=4;
rng(1);
%% Simulation Parameters
M=256;
k=16;
% The number of total antennas
Ms=M/S;
K=k/S; % The number of users
D_vec=2:2:60; % The number of active antennas per user
rho_dB=10; % SNR in dB
rho=10^(rho_dB/10);
P=1; % Total power
Normalization=1; % Normalization type (Possible values 1: Nromalization 1, 2: Normalization 2)
ens=5; % Independent channel realizations for ensemble averaging
for num=1:2
    Normalization=num;
    for s1=1:S
        for s2=1:S
            %% Main loop (over the number of active antennas D)
            for ii=1:length(D_vec)
                D=D_vec(ii);
                %% Best case
                setDbest=([ones(D,K);zeros(Ms-D,K)]); % Initializing active antenna indices for best case
                for kk=1:K
                    setDbest(:,kk)=circshift(setDbest(:,kk),(kk-1)*D); % Circulant shifts to obtain the active antenna indices for best case
                end
                %% Worst case
                setDworst=[ones(D,K);zeros(Ms-D,K)]; % active antenna indices for worst case
                %% Loop for ensemble averaging
                for jj=1:ens               
                    %% Simulation Stationary
                    H=1/sqrt(4)*(randn(Ms,K,s1,s2)+j*randn(Ms,K,s1,s2)); % The channel matrix for stationary case
                    if(Normalization==1)
                        Hbest(:,:,s1,s2)=H(:,:,s1,s2).*setDbest*sqrt(Ms/D); % The channel matrix for the best case
                        Hworst(:,:,s1,s2)=H(:,:,s1,s2).*setDworst*sqrt(Ms/D); % The channel matrix for the worst case
                    elseif(Normalization==2)
                        Hbest(:,:,s1,s2)=H(:,:,s1,s2).*setDbest; % The channel matrix for the best case
                        Hworst(:,:,s1,s2)=H(:,:,s1,s2).*setDworst; % The channel matrix for the worst case
                    end
                end
            end
        end
    end
    if (num==1)
        H_best_nor1=Hbest;
        H_worst_nor1=Hworst;
    else
        H_best_nor2=Hbest;
        H_worst_nor2=Hworst;
    end
end
save('h1_save',"H_worst_nor1","H_best_nor2","H_best_nor1","H_worst_nor2")
