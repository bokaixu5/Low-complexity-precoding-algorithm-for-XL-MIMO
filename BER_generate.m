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
    par.channel = 'rayleigh';   % channel model 'rayleigh'
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
        case 'rayleigh'
             load("h_save.mat","H_best_nor1","H_best_nor2","H_worst_nor1","H_worst_nor2")
            a=0.1;
            tau=0.6; 
            Phi=a.^toeplitz([0:par.B-1]);
            sqrtmPhi = sqrtm(Phi);
            %H = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B));
            W=(randn(par.B,par.U)+1i*randn(par.B,par.U))/sqrt(2*par.U);
    Hn=sqrtmPhi*H_best_nor2;
    error_channel=(randn(par.B,par.U)+1i*randn(par.B,par.U))/sqrt(2*par.U);
    H = (sqrt(1-tau^2)*Hn + tau*sqrtmPhi*error_channel)';
        otherwise
            load("h_save.mat","H_best_nor1","H_best_nor2","H_worst_nor1","H_worst_nor2")
            H=H_best_nor2';
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
    V_RKA = functionRKA(par.B,par.U,1,numRealizations,numIterations(1,b,2),H',updateSchedule(2));
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
numIterations = meanConv;
%Go through all the bounds
for b = 1:2
    %Run RKA
    V_RKA = functionRKA(par.B,par.U,1,numRealizations,numIterations(1,b,2),H',updateSchedule(1));
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
%% Zero-forcing (ZF) beamforming
function [X, beta] = MMSE(par,S,H,N0)

% transmitted signal
P = H'*inv(H*H'+N0/par.rho2*eye(par.U));
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);

% average scaling over signals
beta = 1/betainv;

end

% ZF pseudo inverse
function Hinv = zfinv(par,H)
[U,S] = size(H);
if S>=U
    Hinv = H'/(H*H');
else
    Hinv = (H'*H)\H';
end
end


%% Centralized Wiener-filter (WF) beamforming
function [X, beta] = WF(par,S,H,N0)

% compute regularized inverse
Ainv = rpinv(par,H,par.U*N0/par.rho2);

% calculate precoding factor efficiently
beta = sqrt((par.Es/par.rho2)*(trace(Ainv)-sum(abs(Ainv(:)).^2)*(par.U*N0/par.rho2)));

% apply inverse to data in centralized manner
if par.B>=par.U
    X = (1/beta)*(H'*(Ainv*S));
else
    X = (1/beta)*(Ainv*(H'*S));
end
end
% Wiener filter (WF) regularized pseudo inverse
function Ainv = rpinv(par,H,reg)
[U,S] = size(H);
if S>=U
    Ainv = inv(H*H'+(reg)*eye(par.U));
else
    Ainv = inv(H'*H+(reg)*eye(par.S));
end
end
%% Partially-Decentralized Wiener-filter (WF) beamforming
function [X, beta] = PD_WF(par,S,H,N0)
% decentralized gram matrix computation and localized averaging
Gc = zeros(par.U,par.U);
for cc=1:par.C
    % calculate local Gram matrix
    Hc(:,:,cc) = H(:,(cc-1)*par.S+1:cc*par.S);
    % average among clusters (can be done in a tree-like fashion)
    Gc = Gc + Hc(:,:,cc)*Hc(:,:,cc)';
end
% compute whitening filter at centralized node
Ainv = inv(Gc+(par.U*N0/par.rho2)*eye(par.U));
% calculate precoding factor efficiently
beta = sqrt((par.Es/par.rho2)*(trace(Ainv)-sum(abs(Ainv(:)).^2)*(par.U*N0/par.rho2)));
% whiten transmit signals
Z = (1/beta)*(Ainv*S);
% perform decentralized MRT with whitened signals
for cc=1:par.C
    X((cc-1)*par.S+1:cc*par.S,:) = Hc(:,:,cc)'*Z;
end

end
%% Fully-Decentralized Wiener-filter (WF) beamforming
function [X, beta] = FD_WF(par,S,H,N0)

%initialization
stomp = par.FD_WF.stomp;

% perform fully decentralized WF precoding
for cc=1:par.C

    % calculate local precoding matrix
    Hc = H(:,(cc-1)*par.S+1:cc*par.S);
    Ainvc = rpinv(par,Hc,stomp*par.U*N0/par.rho2);
    betac(cc,1) = sqrt((par.Es/(par.rho2/par.C))*(trace(Ainvc)-sum(abs(Ainvc(:)).^2)*(stomp*par.U*N0/(par.rho2))));

    %betainv(cc,1) = sqrt(par.rho2/par.C)/sqrt(par.Es*trace(Pc*Pc'));

    if par.S>=par.U
        X((cc-1)*par.S+1:cc*par.S,:) = (1/betac(cc,1))*(Hc'*(Ainvc*S));
    else
        X((cc-1)*par.S+1:cc*par.S,:) = (1/betac(cc,1))*(Ainvc*(Hc'*S));
    end

end
% exact beta is hard to compute; just do nothing smart and throw some NaNs
beta = NaN;

end
%% decentralized precoder, ADMM version taken from our old journal paper
function [X,beta] = DP_legacy(par,S,H,N0)

%  -- initialize
delta = par.DP_legacy.delta;
gamma = par.DP_legacy.gamma;
maxiter = par.DP_legacy.maxiter;

H_c = zeros(par.U,par.S,par.C);
AinvH = zeros(par.S,par.U,par.C);

% -- preprocessing
for c=1:par.C
    H_c(:,:,c) = H(:,par.S*(c-1)+1:par.S*c); % get the appropriate part of H
    AinvH(:,:,c) = (H_c(:,:,c)'*H_c(:,:,c) + (1/delta)*(par.U*N0/par.rho2)*eye(par.S))\H_c(:,:,c)';
end
for tt=1:par.T
    % initialize running variables
    lambda_c = zeros(par.U,par.C);
    x_c = zeros(par.S,par.C);
    w_c = zeros(par.U,par.C);
    Hx_c = zeros(par.U,par.C);
    s = S(:,tt);
    % important for fast convergence (reasonable initial guess)
    z_c = max(par.U/par.B,1/par.C)*s*ones(1,par.C);
    % -- start iteration
    for ll = 1:maxiter
        % cluster-wise equalization
        for c=1:par.C
            x_c(:,c) = AinvH(:,:,c)*(z_c(:,c) + lambda_c(:,c));
            Hx_c(:,c) = H_c(:,:,c)*x_c(:,c);
            w_c(:,c) = (Hx_c(:,c)-lambda_c(:,c));
        end
        % consensus step
        w_avg = 1/(par.C*delta+delta^2)*(par.C*s+delta*sum(w_c,2));
        % cluster-wise update
        for c=1:par.C
            z_c(:,c) = (1/delta)*(s+delta*w_c(:,c))-w_avg;
            lambda_c(:,c) = lambda_c(:,c) - gamma*(Hx_c(:,c)-z_c(:,c));
        end
    end
    X(:,tt) = x_c(:); % vectorize output
end
% instantaneous power normalization
X = X*(sqrt(par.rho2)/sqrt(sum(abs(X(:)).^2)/par.T));
% exact beta is hard to compute; just do nothing smart and throw some NaNs
beta = NaN;
end