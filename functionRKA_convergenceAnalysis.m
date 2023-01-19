function conv = functionRKA_convergenceAnalysis(M,K,p,numRealizations,R,H,bounds,updateSchedule)
% Finds the number of iterations necessary to rKA-based schemes converge to a certain performance level.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of antennas.
% @param  K               number of users.
% @param  p               total uplink transmit power per UE [mW].
% @param  numRealizations number of channel realizations (small-fading).
% @param  R               M x M x K matrix with K x (M x M) covariance matrices of K users.
% @param  H               M x numRealizations x K matrix with channel responses.
% @param  bounds          vector with desired performance bounds.
% @param  updateSchedule  different ways of randomly select the equations in rKA. Possible ways: 'power', 'uniform', and 'aa'.
% @return conv            length(bounds) x 1 vector with number of iterations necessary to convergence.
%

%% Preamble

%Compute the inverse of the SNR
xi = 1/p;

%Size of searching region
searchingSize = 1e5;

%Compute RZF sum SINR
sumSINR_RZF = functionComputeSINR(M,K,p,numRealizations,H);

%Prepare to save required number of iterations until convergence
conv = zeros(length(bounds),1);

%% Defining different update schedules

%Prepare to save the update schedule vector
us = zeros(searchingSize,numRealizations);

%Go through all channel realizations
for n = 1:numRealizations

    %Extract channel realizations from all users to BS
    Hn = reshape(H(:,n,:),[M K]);

    %Check update schedule
    if strcmp(updateSchedule,'power') %power-based

        %Define the function to compute the sample probability (SP) based on
        %power ratio
        pwrBasedSP = @(col) ((norm(Hn(:,col))^2 + xi)/((norm(Hn,'fro')^2) + xi*K));

    elseif strcmp(updateSchedule,'uniform') %uniform-based

        %Prepare to save indexes of active users
        activeUsersIndexes = zeros(K,1);

        %Go through all users
        for k = 1:K

            %Compute signal power
            signal = abs(Hn(:,k)'*Hn(:,k)).^2;

            %Check signal power
            if gt(signal, 0)

                %Save index of active user
                activeUsersIndexes(k) = k;

            end

        end

        %Uniform distribution of active users
        unifBasedSP = zeros(K,1);
        unifBasedSP(activeUsersIndexes > 0) = 1/sum(activeUsersIndexes > 0);

    elseif strcmp(updateSchedule,'aa')

        %Prepare to save number of active antennas of each user
        activeAntennas = zeros(K,1);

        %Go through all users
        for k = 1:K

            %Compute number of active antennas of user k
            activeAntennas(k) = trace(R(:,:,k) > 0);

        end

        %Define sample probabilities based on active antennas
        aaBasedSP = activeAntennas/sum(activeAntennas);

    end

    %Check and generate the update schedule (us)
    if strcmp(updateSchedule,'power')

        us(:,n) = functionRandp(arrayfun(pwrBasedSP,(1:K)'),searchingSize,1);

    elseif strcmp(updateSchedule,'uniform')

        us(:,n) = functionRandp(unifBasedSP,searchingSize,1);

    elseif strcmp(updateSchedule,'aa')

        us(:,n) = functionRandp(aaBasedSP,searchingSize,1);

    end

end

%% Obtaining the average number of iterations until convergence

%Start number of iterations
int = 1;

%Start bound index
bindex = 1;

%Define stopping vector
stopVector = zeros(length(bounds),1);

%Define looping stop criteria
stop = false;

%Prepare to store state vectors
u = zeros(M,2,numRealizations,K);
z = zeros(K,2,numRealizations,K);

%Prepare to save receive combining matrix
V = zeros(M,numRealizations,K);

%Go through all iterates
while not(stop)

    %Go through all channel realizations
    for n = 1:numRealizations

        %Extract channel realizations from all users to BS
        Hn = reshape(H(:,n,:),[M K]);

        %Go through all users
        for k = 1:K

            %Compute signal energy
            signal = p*abs(Hn(:,k)'*Hn(:,k)).^2;

            %Check if the energy of the signal is nonzero
            if gt(signal,0)

                %Define canonical basis
                e = zeros(K,1); e(k) = 1;

                if int == 1 %self-initialization

                    eqPick = k;

                else

                    eqPick = us(int,n);

                end

                %Pick channel vector
                q = Hn(:,eqPick);

                %Compute residual
                eta = (e(eqPick) - (q'*u(:,1,n,k)) - xi*z(eqPick,1,n,k))/((norm(q)^2) + xi);

                %Updating state vectors
                u(:,2,n,k) = u(:,1,n,k) + eta*q;
                z(:,2,n,k) = z(:,1,n,k);
                z(eqPick,2,n,k) = z(eqPick,1,n,k) + eta;

                %Compute receive combining matrix
                V(:,n,k) = Hn*z(:,1,n,k);

            end

        end

    end

    %Compute RKA sum SINR
    sumSINR_RKA = functionComputeSINR(M,K,p,numRealizations,H,V);

    if ge((sumSINR_RKA/sumSINR_RZF),bounds(bindex))

      stopVector(bindex) = 1;
      conv(bindex) = int;
      bindex = bindex + 1;

    end

    if sum(stopVector) == length(bounds)

        stop = true;

    else %continuing

        %Saving the context
        u(:,1,:,:) = u(:,2,:,:);
        z(:,1,:,:) = z(:,2,:,:);

        %Incrementing the number of iterations
        int = int + 1;

    end

end

end
