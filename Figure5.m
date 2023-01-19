close all;
clear;
%Select length of coherence block
tau_c = 250;
%Define range of number of BS antennas
Mrange = 10:1:256;
%Define range of number of UEs
Krange = 1:100;
%Set number of cells considered in the M-MMSE scheme
L = 9;
numIter=100;
%% Consider M=100 and varying K
K = Krange;
M = max(Mrange);
%Compute number of samples for uplink data
tau_u = (tau_c-K);
%Compute complexity of receive combining
receiverProcessing = tau_u.*K*M;
%Compute rKA computational complexity
rKA=receiverProcessing+M*numIter + M;
s_rKA=receiverProcessing+M*numIter + 2*M*K;
%Add complexity of computing combining matrix
complexity_MMMSE = receiverProcessing + L*K*(M^2+M)/2 + M^2*K + (M^3-M)/3; %M-MMSE
complexity_SMMSE = receiverProcessing + 3*M^2*K/2 + M*K/2 + (M^3-M)/3; %S-MMSE
complexity_RZF = receiverProcessing + 3*K.^2*M/2 + 3*M*K/2 + (K.^3-K)/3; %S-MMSE2
complexity_ZF = receiverProcessing + 3*K.^2*M/2 + M*K/2 + (K.^3-K)/3; %ZF
complexity_MR = receiverProcessing; %MR
%Plot the simulation results for M=100 and varying K
figure(1);
hold on; box on;
%plot(K(1),complexity_MMMSE(1),'rd-','LineWidth',1);
%plot(K(1),complexity_SMMSE(1),'b:','LineWidth',2);
plot(K(1),complexity_RZF(1),'k-.','LineWidth',3);
%plot(K(1),complexity_ZF(1),'r--','LineWidth',1);
%plot(K(1),complexity_MR(1),'bs-','LineWidth',1);
plot(K,rKA,'b-','LineWidth',3);
plot(K,s_rKA,'r--','LineWidth',3);
%plot(K,complexity_MMMSE,'r-','LineWidth',1);
%plot(K,complexity_SMMSE,'b:','LineWidth',2);
plot(K,complexity_RZF,'k-.','LineWidth',3);
%plot(K,complexity_ZF,'r--','LineWidth',1);
plot(K,rKA,'b-','LineWidth',3);
plot(K,s_rKA,'r--','LineWidth',3);
%plot(K,complexity_MR,'b-','LineWidth',1);
%plot(K([1 5:5:40]),complexity_MMMSE([1 5:5:40]),'rd','LineWidth',1);
%plot(K([1 5:5:40]),complexity_MR([1 5:5:40]),'bs','LineWidth',1);
xlabel('Number of UEs (K)','Interpreter','Latex');
ylabel('Number of complex multiplications','Interpreter','Latex');
set(gca,'YScale','log');
set(gca,'FontSize',12);
set(gca,'xLim',[50 100]);
legend('RZF','rKA','SwoR-rKA','Location','SouthEast','Interpreter','Latex');
grid on
%% Consider K=10 and varying M
K = Krange(10);
M = Mrange;
%Compute number of samples for uplink data
tau_u = (tau_c-K);
%Compute complexity of receive combining
receiverProcessing = tau_u*K*M;
%Add complexity of computing combining matrix
complexity_MMMSE = receiverProcessing + L*K*(M.^2+M)/2 + M.^2*K + (M.^3-M)/3; %M-MMSE
complexity_SMMSE = receiverProcessing + 3*M.^2*K/2 + M*K/2 + (M.^3-M)/3; %S-MMSE
complexity_RZF = receiverProcessing + 3*K.^2*M/2 + 3*M*K/2 + (K^3-K)/3; %RZF
complexity_ZF = receiverProcessing + 3*K^2*M/2 + M*K/2 + (K^3-K)/3; %ZF
complexity_MR = receiverProcessing; %MR
rKA=receiverProcessing+M*numIter + M;
s_rKA=receiverProcessing+M*numIter + 2*M*K;
%Plot the simulation results for K=10 and varying M
figure(2);
hold on; box on;
%plot(M(1),complexity_MMMSE(1),'rd-','LineWidth',1);
%plot(M(1),complexity_SMMSE(1),'b:','LineWidth',3);
plot(M(1),complexity_RZF(1),'k-.','LineWidth',3);
%plot(M(1),complexity_ZF(1),'r--','LineWidth',1);
%plot(M(1),complexity_MR(1),'bs-','LineWidth',1);
plot(M,rKA,'b-','LineWidth',3);
plot(M,s_rKA,'r--','LineWidth',3);
%plot(M,complexity_MMMSE,'r-','LineWidth',1,);
%plot(M,complexity_SMMSE,'b:','LineWidth',2);
plot(M,complexity_RZF,'k-.','LineWidth',3);
%plot(M,complexity_ZF,'r--','LineWidth',1);
%plot(M,complexity_MR,'b-','LineWidth',1);
plot(M,rKA,'b-','LineWidth',3);
plot(M,s_rKA,'r--','LineWidth',3);
%plot(M,complexity_MMMSE,'rd','LineWidth',1);
%plot(M,complexity_MR,'bs','LineWidth',1);
xlabel('Number of antennas (M)','Interpreter','Latex');
ylabel('Number of complex multiplications','Interpreter','Latex');
set(gca,'YScale','log');
%set(gca,'YLim',[50000 5000000]);
set(gca,'xLim',[100 200]);
legend('RZF','rKA','SwoR-rKA','Location','NorthWest','Interpreter','Latex');
set(gca,'FontSize',12);
grid on