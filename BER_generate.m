%Initialization
close all;
clearvars;
% -- set up default/custom parameters
    disp('using default simulation settings and parameters...')
    % set default simulation parameters
    par.runId = 0;              % simulation ID (used to reproduce results)
    par.U = 16;                 % number of single-antenna users
    par.B = 256;                % number of base-station antennas (B>>U)
    par.T = 10;                 % number of time slots
    par.C = 16;                 % number of antenna clusters
    par.mod = '64QAM';          % modulation type: 'BPSK','QPSK','16QAM','64QAM','8PSK'
    par.trials = 1e2;           % number of Monte-Carlo trials (transmissions)
    par.NTPdB_list = -10:2:40;  % list of normalized transmit power [dB] values
    par.rho2 = 1;               % rho^2=1 (should NOT affect your results!)
    par.precoder = {'RZF','rKA','SwoR-rKA'};
    par.channel = 'Imperfect CSI XL-MIMO Channel';   % channel model 'Imperfect CSI XL-MIMO Channel'
    par.betaest = 'pilot';      % 'pilot'
    par.save = false;            % save results (true,false)
    par.plot = true;            % plot results (true,false)
    % algorithm-dependent parameters
    % FD_WF (please tune if you change the scenario)
    par.FD_WF.stomp = 0.125*par.C; % determines how much to stomp regularization (C)
    % DP_legacy (please tune if you change the scenario)
    par.DP_legacy.delta = 0.3; % Lagrange scaling (1.0)
    par.DP_legacy.gamma = 1.0; % Lagrange stepsize (1.0)
    par.DP_legacy.maxiter = 2; % keep at 2
% -- initialization
% use runId random seed (enables reproducibility)
rng(par.runId);
% simulation name (used for saving results)
par.simName = ['ERR_',num2str(par.U),'x',num2str(par.B), '_C', ...
    num2str(par.C), '_', par.betaest, '_', par.mod, '_', num2str(par.trials),'Trials'];
% set up Gray-mapped constellation alphabet (some are selected according to IEEE 802.11)
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    case '8PSK',
        par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
            exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
            exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
            exp(1i*2*pi/8*4), exp(1i*2*pi/8*5) ];
    case '16PSK',
        par.symbols = [ ...
            exp(1i*2*pi*0/16)  ... % 0000
            exp(1i*2*pi*1/16)  ... % 0001
            exp(1i*2*pi*3/16)  ... % 0010
            exp(1i*2*pi*2/16)  ... % 0011
            exp(1i*2*pi*7/16)  ... % 0100
            exp(1i*2*pi*6/16)  ... % 0101
            exp(1i*2*pi*4/16)  ... % 0110
            exp(1i*2*pi*5/16)  ... % 0111
            exp(1i*2*pi*15/16) ... % 1000
            exp(1i*2*pi*14/16) ... % 1001
            exp(1i*2*pi*12/16) ... % 1010
            exp(1i*2*pi*13/16) ... % 1011
            exp(1i*2*pi*8/16)  ... % 1100
            exp(1i*2*pi*9/16)  ... % 1101
            exp(1i*2*pi*11/16) ... % 1110
            exp(1i*2*pi*10/16) ];  % 1111
end

% compute symbol energy
par.Es = mean(abs(par.symbols).^2);

% number of antenns per cluster
par.S = par.B/par.C;

% precompute bit labels
par.bps = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.bps,'left-msb');

% track simulation time
time_elapsed = 0;
j=sqrt(-1);
rng(1);
%% Simulation Parameters
M=256; % The number of total antennas
K=16; % The number of users
D_vec=2:2:60; % The number of active antennas per user
rho_dB=10; % SNR in dB
rho=10^(rho_dB/10);
P=1; % Total power
Normalization=2; % Normalization type (Possible values 1: Nromalization 1, 2: Normalization 2)
%
ens=1e2; % Independent channel realizations for ensemble averaging
%% Main loop (over the number of active antennas D)
for ii=1:length(D_vec)
    D=D_vec(ii);
    %% Best case
    setDbest=[ones(D,K);zeros(M-D,K)]; % Initializing active antenna indices for best case
    for kk=1:K
        setDbest(:,kk)=circshift(setDbest(:,kk),(kk-1)*D); % Circulant shifts to obtain the active antenna indices for best case
    end
    %% Worst case
    setDworst=[ones(D,K);zeros(M-D,K)]; % active antenna indices for worst case
    %% Loop for ensemble averaging   
        %% Simulation Stationary
        H=1/sqrt(2)*(randn(M,K)+j*randn(M,K)); % The channel matrix for stationary case
        if(Normalization==1)
            Hbest=H.*setDbest*sqrt(M/D); % The channel matrix for the best case
            Hworst=H.*setDworst*sqrt(M/D); % The channel matrix for the worst case
        elseif(Normalization==2)
            Hbest=H.*setDbest; % The channel matrix for the best case
            Hworst=H.*setDworst; % The channel matrix for the worst case
        end
end    
% -- start simulation

% - initialize result arrays (detector x normalized transmit power)
[res.PER, res.SER, res.BER ] = deal(zeros(length(par.precoder),length(par.NTPdB_list)));
[res.TxPower, res.RxPower, res.TIME] = deal(zeros(length(par.precoder),length(par.NTPdB_list)));

% compute noise variances to be considered: NTP = rho^2/N0
N0_list = par.rho2*10.^(-par.NTPdB_list/10);
% trials loop
tic
for t=1:par.trials

    % generate data
    for qq=1:par.T
        % generate random bit stream
        B(:,:,qq) = randi([0 1],par.U,par.bps);
        % generate transmit symbol
        Idx(:,qq) = bi2de(B(:,:,qq),'left-msb')+1;
        S(:,qq) = par.symbols(Idx(:,qq)).';
    end

    % generate masks
    MaskI = true(par.U,par.T);
    MaskT = true(1,par.T);
    % here you could add other estimation methods
    switch par.betaest
        case {'pilot'}
            S(:,1) = ones(par.U,1)*sqrt(par.Es); % just send all ones in the first time slots
            MaskI(:,1) = false;
            MaskT(1,1) = false;
        otherwise,
            error('par.betaest not specified')
    end

    % generate iid Gaussian channel matrix and noise matrix
    N = sqrt(0.5)*(randn(par.U,par.T)+1i*randn(par.U,par.T));

    % you can add your own channel model here
    switch par.channel
        case 'Imperfect CSI XL-MIMO Channel'
            a=0.1;
            tau=0.3; 
            Phi=a.^toeplitz([0:par.B-1]);
            sqrtmPhi = sqrtm(Phi);
            %H = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B));
            W=(randn(par.B,par.U)+1i*randn(par.B,par.U))/sqrt(2*par.U);
    Hn=sqrtmPhi*Hbest;
    error_channel=(randn(par.B,par.U)+1i*randn(par.B,par.U))/sqrt(2*par.U);
    H = (sqrt(1-tau^2)*Hn + tau*sqrtmPhi*error_channel)';
        otherwise
            %Perfect CSI XL-MIMO Channel
            H=Hbest';
    end

    % algorithm loop
    for d=1:length(par.precoder)

        % normalized transmit power loop
        for k=1:length(par.NTPdB_list)

            % set noise variance
            N0 = N0_list(k);

            % record time used by the beamformer
            starttime = toc;

            % beamformers
            switch (par.precoder{d})
               
                case 'rKA',
                    [X,beta]=RKA(par,S,H,N0);
                case 'SwoR-rKA',
                    [X,beta]=RKA2(par,S,H,N0);
                case 'RZF',
                    [X,beta]=RZF(par,S,H,N0);
                    case 'ZF',
                    [X,beta]=ZF(par,S,H,N0);
               
                otherwise,
                    error('par.precoder not specified')

            end

            % record beamforming simulation time
            res.TIME(d,k) = res.TIME(d,k) + (toc-starttime);


            % transmit data over noisy channel
            HX = H*X;
            Y = HX + sqrt(N0)*N;

            % extract transmit and receive power
            res.TxPower(d,k) = res.TxPower(d,k) + mean(sum(abs(X(:)).^2))/par.T;
            res.RxPower(d,k) = res.RxPower(d,k) + mean(sum(abs(HX(:)).^2))/par.U/par.T;

            % UEs must estimate beta
            switch par.betaest
                case 'genie', % perfect beta directly from beamformer
                    betaest = ones(par.U,1)*beta;
                case 'pilot', % knows that first symbols are for training
                    betaest = real(1./Y(:,1)*sqrt(par.Es)); % ML estimate since we have no prior on beta
            end

            % perform estimation
            Shat = (betaest*ones(1,par.T)).*Y;

            % UE-side hard-output data detection
            for qq=1:par.T
                [~,Idxhat(:,qq)] = min(abs(Shat(:,qq)*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
                Bhat(:,:,qq) = par.bits(Idxhat(:,qq),:);
            end

            % -- compute error and complexity metrics
            err = (Idx(MaskI)~=Idxhat(MaskI));
            res.PER(d,k) = res.PER(d,k) + any(err(:));
            res.SER(d,k) = res.SER(d,k) + sum(err(:))/par.U/par.T;
            tmpBER = B(:,:,MaskT)~=Bhat(:,:,MaskT);
            res.BER(d,k) = res.BER(d,k) + sum(tmpBER(:))/(par.U*par.bps*sum(MaskT));

        end % NTP loop

    end % algorithm loop

    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',...
            time_elapsed*(par.trials/t-1)/60);
        tic
    end

end % trials loop

% normalize results
res.PER = res.PER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.TxPower = res.TxPower/par.trials;
res.RxPower = res.RxPower/par.trials;
res.TIME = res.TIME/par.trials;

% manually (or visually) check whether the TX power of your precoder is
% correct (this is a very common mistake in many papers...)
res.TxPower

% -- save final results (par and res structures)

if par.save
    save([ par.simName '_' num2str(par.runId) ],'par','res');
end

% -- show results (generates fairly nice Matlab plots)

if par.plot

    % - BER results
    marker_style = {'kx-','bo:','rs--','mv-.','gp-.','bs--','y*--'};
    h = figure(1);
    for d=1:length(par.precoder)
        semilogy(par.NTPdB_list,res.BER(d,:),marker_style{d},'LineWidth',2);
        if (d==1)
            hold on
        end
    end
    hold off
    grid on
    box on
    xlabel('Normalized transmit power [dB]','FontSize',12,'Interpreter','Latex')
    ylabel('Uncoded bit error rate (BER)','FontSize',12,'Interpreter','Latex');
    if length(par.NTPdB_list) > 1
        axis([min(par.NTPdB_list) max(par.NTPdB_list) 1e-3 1]);
    end
    legend(par.precoder,'FontSize',12,'Interpreter','Latex','location','northeast')
    set(gca,'FontSize',12);
    if par.save
        % save eps figure (in color and with a reasonable bounding box)
        print(h,'-loose','-depsc',[ par.simName '_' num2str(par.runId) ])
    end
end
save('ber3',"par","res")

function [X, beta] = RKA(par,S,H,N0)
numIterations = 100;
updateSchedule = ["power","uniform","aa"];
numRealizations = 20;
%Go through all the bounds
for b = 1:2
    %Run RKA
    V_RKA = functionRKA(par.B,par.U,1,numRealizations,numIterations,H',updateSchedule(2));
end
P=(reshape(V_RKA(:,3,:),[par.B par.U]));
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);
% average scaling over signals
beta = 1/betainv;
end

function [X, beta] = RKA2(par,S,H,N0)
numIterations = 100;
updateSchedule = ["power","uniform","aa"];
numRealizations = 20;
%Go through all the bounds
for b = 1:2
    %Run RKA
    V_RKA = functionRKA(par.B,par.U,1,numRealizations,numIterations,H',updateSchedule(1));
end
P=(reshape(V_RKA(:,3,:),[par.B par.U]));
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);
% average scaling over signals
beta = 1/betainv;
end
%% Maximum ratio transmission (MRT) beamforming
function [X, beta, P] = MRT(par,S,H,N0)

% transmitted signal
P = H';
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);

% average scaling over signals
beta = 1/betainv;

%beta = 1.0*(norm(s,2)^2+N0*par.U)/(s'*H*x); % that's cheating

end

function [X, beta] = ZF(par,S,H,N0)

% transmitted signal
P = zfinv(par,H);
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);

% average scaling over signals
beta = 1/betainv;

end
function [X, beta] = RZF(par,S,H,N0)

% transmitted signal
P =  H'*inv(H*H'+ par.U/par.rho2*eye(par.U));
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);

% average scaling over signals
beta = 1/betainv;

end
