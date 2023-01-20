%Initialization
close all;
clearvars;
%Fixing random seed
rng(0)
%Number of antennas
M = 256;
S=4;
Ms=M/S;
%Number of users
K=16;
K = K/S;
%Normalization of the diagonal matrix
diagNorm = 'Norm2';
updateSchedule = ["power","uniform","aa"];
%Set number of setups
numSetups = 3;
%Number of stats ergodic realizations
numRealizations = 20;
%% Loading simulation setup
j=sqrt(-1);
rng(1);
%% Simulation Parameters
D_vec=2:2:60; % The number of active antennas per user
rho_dB=10; % SNR in dB
rho=10^(rho_dB/10);
P=1; % Total power
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
%% Algorithm parameters
%Number of bounds
bounds = 2;
Hnum=100;
snr = -10:3:20;
sumratedataRZF = zeros(1,length(snr));
sumratedataRKA=zeros(1,length(snr));
sumratedataRKA2=zeros(1,length(snr));
a=0.1;						% Covariance/Correlation between antennas
tau=0.1;                    % CSI quality
% Monte-Carlo runs
num_iter=10;               % Choose >1000 for print quality
%% Channel Covariance Matrix
Phi=a.^toeplitz([0:M-1]);
sqrtmPhi = sqrtm(Phi);
time_elapsed = 0;
tic;
%Go through all different SNR values
for idx1 = 1:1:length(snr)
    SRrzf = 0;
    SRzf = 0;
    SRmrt = 0;
    SRwmmse = 0;
    SRKA=0;
    p=10^(snr(idx1)/10);
    sigma2=1;
    H=H_best_nor2;
    %Hn = randn(K,Ms,S) + 1i*randn(K,Ms,S);
    %H = sqrt(1/4)*(randn(K,Ms,S,S) + 1i*randn(K,Ms,S,S));
    %% Create Channel draw with error and covariance
    %     W=(randn(M,K)+1i*randn(M,K))/sqrt(2*K);
    %     H=sqrtmPhi*W;
    %H=H_best_nor2;
    %     h_1=H(:,1);
    %     error_channel=(randn(M,K)+1i*randn(M,K))/sqrt(2*K);
    %     H_hat = sqrt(1-tau^2)*H + tau*sqrtmPhi*error_channel;
    %Hn=H_hat';
    %Hn=H_best_nor2';
    %Go through all the bounds
    for idx2 = 1:1:Hnum
        a=1;
        SRrzf = SRrzf + SumRate(H,a,sigma2,S,p);
        a=2;
        SRKA = SRKA + SumRate(H,a,sigma2,S,p);
        a=3;
        SRKA2 = SRKA + SumRate(H,a,sigma2,S,p);
    end
    sumratedataRZF(idx1) = SRrzf/Hnum;
    sumratedataRKA(idx1) = SRKA /Hnum;
    sumratedataRKA2(idx1) = SRKA2 /Hnum;
    time=toc;
    time_elapsed = time_elapsed + time;
    fprintf('estimated remaining simulation time: %3.0f min.\n',...
        time_elapsed*(length(snr)/idx1-1)/60);
    tic
end
figure
plot(snr,sumratedataRZF,'k-*','LineWidth',2)
%plot(snr,sumratedataWMMSE,'b-^','LineWidth',4)
hold on;
plot(snr,sumratedataRKA,'R:^','LineWidth',2)
plot(snr,sumratedataRKA2,'b-.','LineWidth',2)
set(gca,'FontSize',12);
set(gca,'xLim',[5,20]);
xlabel('SNR [dB]','Interpreter','Latex')
ylabel('Average sum SE [bit/s/Hz/subarray]','Interpreter','Latex')
legend('RZF','rKA','SwoR-rKA','Location','best','Interpreter','Latex')
grid on