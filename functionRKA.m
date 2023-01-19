function [V,us] = functionRKA(M,K,p,numRealizations,numInterRange,H,updateSchedule)
% Outputs the receive combining matrix V and the probabilities probUS related to the update schedule of rKA-based schemes.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of antennas.
% @param  K               number of users.
% @param  p               total uplink transmit power per UE [mW].
% @param  numRealizations number of channel realizations (small-fading).
% @param  numIterRange    vector with number of iterations points.
% @param  R               M x M x K matrix with K x (M x M) covariance matrices of K users.
% @param  H               M x numRealizations x K matrix with channel responses.
% @param  updateSchedule  different ways of randomly select the equations in rKA. Possible ways: 'power', 'uniform', and 'aa'.
% @return V               M x numRealizations x K receive combining matrix.
% @return probUS          K x K matrix with users's selection probabilities.
%

%% Preamble

%Compute the inverse of the SNR
xi = 1/p;

%Extract the max number of iterations
maxNumIterations = 200;

%Prepare to store receive combining matrix
V = zeros(M,numRealizations,K,length(numInterRange));

%Prepare to store probabilities of the update schedule
probUS = zeros(K,K);

%% Go through all channel realizations
for n = 1:numRealizations

    %Extract channel realizations from all users to the BS
    Hn = H;

    %Check update schedule
    if strcmp(updateSchedule,'power')

        %Define the function to compute the sample probability (SP) based on
        %the power ratio, a weigth that considers the pwr of the equation in
        %relation to the SLE
        pwrBasedSP = @(col) ((norm(Hn(:,col))^2 + xi)/((norm(Hn,'fro')^2) + xi*K));

    elseif strcmp(updateSchedule,'uniform')

        %Prepare to save the indexes of active users
        activeUsersIndexes = zeros(K,1);

        %Go through all users
        for k = 1:K

            %Compute the signal power
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

        %Prepare to save the quantity of active antennas of each users
        activeAntennas = zeros(K,1);

        %Go through all users
        for k = 1:K

            %Compute the number of active antennas of user k
            activeAntennas(k) = trace(R(:,:,k) > 0);

        end

        %Define sample probabilities based on active antennas
        aaBasedSP = activeAntennas/sum(activeAntennas);

    end

    %Go through all users 
    for k = 1:K

        %Compute signal energy
        signal = p*abs(Hn(:,k)'*Hn(:,k)).^2;

        %Check for nonzero signal energy
        if gt(signal,0)

            %Define canonical basis
            e = zeros(K,1); e(k) = 1;

            %Prepare to store state vectors
            u = zeros(M,maxNumIterations);
            z = zeros(K,maxNumIterations);

            %Check and generate the update schedule (us)
            if strcmp(updateSchedule,'power')

                us = functionRandp(arrayfun(pwrBasedSP,(1:K)'),maxNumIterations,1);

            elseif strcmp(updateSchedule,'uniform')

                us = functionRandp(unifBasedSP,maxNumIterations,1);

            elseif strcmp(updateSchedule,'aa')

                us = functionRandp(aaBasedSP,maxNumIterations,1);

            elseif strcmp(updateSchedule,'sir')

                us = functionRandp(sirBasedSP(:,k),maxNumIterations,1);

            end

            %Go through all users that can be selected
            for kk = 1:K

                %Storing the probability of selection
                probUS(kk,k) = probUS(kk,k) + mean(us == kk)/numRealizations;

            end

            %Go through all iterations values
            for int = 1:maxNumIterations

                if int == 1 %self-initialization

                    eqPick = k;

                else

                    eqPick = us(int);

                end

                %Pick channel vector
                q = Hn(:,eqPick);

                %Compute residual
                eta = (e(eqPick) - (q'*u(:,int)) - xi*z(eqPick,int))/((norm(q)^2) + xi);

                %Update schedule
                u(:,int+1) = u(:,int) + eta*q;
                z(:,int+1) = z(:,int);
                z(eqPick,int+1) = z(eqPick,int+1) + eta;

                %Check if number of iterations reached a desired number
                if sum(int == numInterRange)

                    %Store results given the index value
                    V(:,n,k,int == (numInterRange)) = Hn*z(:,int);

                end

            end

        end

    end

end

end
