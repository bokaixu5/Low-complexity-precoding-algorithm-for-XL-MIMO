function y=fixed_point_solve(Phi,tau,rho_lin,K)
%The definitions and equations refer (mostly) to the 
%following publication:
%
%Müller, A., A. Kammoun, E. Björnson, and M. Debbah, "Linear Precoding 
%Based on Polynomial Expansion: Reducing Complexity in Massive MIMO", 
%EURASIP Journal on Wireless Communications and Networking, 2016(1), 1-22, 
%DOI: 10.1186/s13638-016-0546-z.
%Pre-print available at: http://www.laneas.com/axel-muller
%
%This is version 1.0. (Last edited: 2016-03-16)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

%% Starting Point
csi=0.5;    % any number > 0
y=equation_solve(Phi,csi,tau,rho_lin,K);

%% Fixpoint Iteration
while(abs(y-csi)/csi>1e-4)  % Check for convergence/accuracy; 1e-4 should be enough
    csi=y;
    y=equation_solve(Phi,csi,tau,rho_lin,K);
end


end %fcn