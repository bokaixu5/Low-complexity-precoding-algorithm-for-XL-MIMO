function V = functionRZF(M,K,p,numRealizations,H)
% Computes receive combining matrix V based on the regularized zero-forcing scheme.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of BS antennas.
% @param  K               number of users.
% @param  p               total uplink transmit power per UE [mW].
% @param  numRealizations number of channel realizations (small-fading).
% @param  H               M x numRealizations x K matrix with channel responses.
% @return V               M x numRealizations x K RZF receive combining matrix.
%

%Store the inverse of the uplink transmit power
xi = 1/p;

%Store K x K identity matrix
eyeK = eye(K);

%Prepare to save RZF receive combining matrix
V = zeros(M,numRealizations,K);

%Go through all channel realizations
for n = 1:numRealizations

    %Extract current channel matrix
    Hn = reshape(H(:,n,:),[M K]);

    %Compute RZF receive combining matrix
    V(:,n,:) = reshape(Hn/(Hn'*Hn + xi*eyeK),[M 1 K]);

end
