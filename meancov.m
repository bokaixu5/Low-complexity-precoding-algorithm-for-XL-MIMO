% Matlab script used to generate Figure 2a of the article:
%
% "Low-Complexity Distributed XL-MIMO for Multiuser Detection"
%
% @author Victor Croisfelt Rodrigues
% @author Abolfazl Amiri
% @author Taufik Abrao
% @author Elisabeth de Carvalho
% @author Petar Popovski
%

%Initialization
close all;
clearvars;

%Fixing random seed
rng(0)

%% Define communication setup

%Number of antennas
M = 128;

%Number of subarrays
S = 4;

%Number of antennas per subarray
Ms = M/S;

%Number of users
K = 16;

%Normalization of the diagonal matrix
diagNorm = 'Norm2';

%% Propagation parameters

%Define uplink transmit power [mW]
p = 1;

%Define receiver noise power [dBm]
noiseVariancedBmRange =-50;

%% Algorithm parameters

%Define update schedule vector
updateSchedule = ["power","uniform","aa"];

%Define performance bounds
bounds = [.90,.99];
%bounds = .9;

% IMPORTANT
% Feel free to try different inputs for bounds and updateSchedule (considering "power","uniform","aa"). 
% However, the plot code in the end of this file is expecting the use of: 
%     bounds = [.90,.99];
%     updateSchedule = ["power","uniform","aa"];
%
% If you use any other input, you will need to correct the code accordingly.

%% Simulation parameters

%Set number of setups
numSetups = 3;

%Number of stats ergodic realizations
numRealizations = 20;

%% Simulation

%Prepare to save simulation results
meanUsers = zeros(numSetups,1); % average number of users being served per subarray
conv = zeros(length(noiseVariancedBmRange),length(bounds),length(updateSchedule),S,numSetups); % iterations till (performance) convergence

%Go through all setups
for s = 1:numSetups

    setup = tic; % setting up a timer: setup
    warning('off','MATLAB:sqrtm:SingularMatrix') % omitting warning related to sqrtm function

    %Output setup progress
    disp([num2str(s) ' setup out of ' num2str(numSetups)]);

    %Prepare to save parfor results
    par_conv = zeros(length(noiseVariancedBmRange),length(bounds),length(updateSchedule),S);

    %Generate mobile communication setup
    [channelGaindB,R,par_meanUsers] = functionExampleSetup(M,S,K,diagNorm);

    %Go through all different SNR values
    for r = 1:length(noiseVariancedBmRange)

        %Output SNR progress
        disp([num2str(r) ' SNR out of ' num2str(length(noiseVariancedBmRange))]);

        %Extract current noise variance
        noiseVariancedBm = noiseVariancedBmRange(r);

        %Compute channel gain over noise
        channelGainOverNoise = zeros(Ms,K,S);
        channelGainOverNoise(channelGaindB ~= 0) = channelGaindB(channelGaindB ~= 0) - noiseVariancedBm;

        %Go through each subarray
        for ss = 1:S

            %Channel realizations
            H = functionChannelRealizations(Ms,K,channelGainOverNoise(:,:,ss),R(:,:,:,ss),numRealizations);

            %Go through all different USs
            for us = 1:length(updateSchedule)

                %Run randomized Kaczmarz based detection
                par_conv(r,:,us,ss) = functionRKA_convergenceAnalysis(Ms,K,p,numRealizations,R(:,:,:,ss),H,bounds,updateSchedule(us));

            end

        end

    end

    %Save parfor results
    meanUsers(s) = par_meanUsers;
    conv(:,:,:,:,s) = par_conv;

    toc(setup)

end

%% Upper bounds

%Compute average number of users per subarray
Kbar = ceil(mean(meanUsers));

%Compute upper bounds
TupperBoundPower = ceil((1/3)*(Kbar^3)/Ms + (2/3)*Kbar/Ms + (3/2)*Kbar^2 - (1/2)*Kbar);
TupperBoundUniform = ceil((1/3)*(Kbar^3)/Ms + (2/3)*Kbar/Ms + (3/2)*Kbar^2 + (3/2)*Kbar - 1);

%% Data extraction
meanConv = ceil(mean(mean(conv,5),4));

%% Computational Relaxation Degree (CRD)
CRD(:,:,1) = functionCRD(meanConv(:,:,1),TupperBoundPower);
CRD(:,:,2) = functionCRD(meanConv(:,:,2),TupperBoundUniform);
CRD(:,:,3) = functionCRD(meanConv(:,:,3),TupperBoundUniform);

%% Save simulation results
save([date,'_CRD_',diagNorm,'_M',num2str(M),'_S',num2str(S),'_K',num2str(K)])

%% Plotting simulation results

figure(1);
subplot(1,2,2)
hold on; box on; ax = gca;

plot(noiseVariancedBmRange,zeros(length(noiseVariancedBmRange),1),'k-','LineWidth',2)

plot(noiseVariancedBmRange,CRD(:,1,1),'k*','LineWidth',1.5)
plot(noiseVariancedBmRange,CRD(:,1,2),'kd','LineWidth',1.5)
plot(noiseVariancedBmRange,CRD(:,1,3),'kh','LineWidth',1.5)

ax.ColorOrderIndex = 1;

plot(noiseVariancedBmRange,CRD(:,1,1),'*--','LineWidth',1.5)
plot(noiseVariancedBmRange,CRD(:,1,2),'d--','LineWidth',1.5)
plot(noiseVariancedBmRange,CRD(:,1,3),'h--','LineWidth',1.5)

ax.ColorOrderIndex = 1;

plot(noiseVariancedBmRange,CRD(:,2,1),'*-.','LineWidth',1.5)
plot(noiseVariancedBmRange,CRD(:,2,2),'d-.','LineWidth',1.5)
plot(noiseVariancedBmRange,CRD(:,2,3),'h-.','LineWidth',1.5)

xlabel('Noise variance [dBm]');
ylabel('Average computational relaxation degree per subarray, $\overline{\mathrm{CRD}}$');

title('Normalization 1 of $\mathbf{D}_{k}$')

legend('RZF benchmark','Power','Uniform','Active-antennas','Location','best');
