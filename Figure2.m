%Initialization
close all;
clearvars;
%Fixing random seed
rng(0)
%Number of antennas
M = 256;
S=4;
Ms=64;
%Number of users
K = 16;
%Normalization of the diagonal matrix
diagNorm = 'Norm2';
updateSchedule = ["power","uniform","aa"];
%Set number of setups
numSetups = 3;
%Number of stats ergodic realizations
numRealizations = 20;
%% Loading simulation setup
load("h1_save.mat","H_best_nor1","H_best_nor2","H_worst_nor1","H_worst_nor2")
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