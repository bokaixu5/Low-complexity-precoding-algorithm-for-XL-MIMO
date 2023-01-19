function y=equation_solve(Phi,csi,tau,rho_lin,K)
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

M=length(Phi);
t=1/csi;
delta=deterministic(Phi,t,K);

T=inv(eye(M)+t/(1+t*delta)*Phi);
gamma=1/K*trace(Phi^2*T^2);
delta_bar=1/csi*delta;
T_bar=1/csi*T;

mu=csi*1/K*trace(Phi*T^3)/(gamma*1/K*trace(Phi*T^2))*(gamma/(1/K*trace(Phi*T^2))-1/K*trace(Phi^2*T^3)/(1/K*trace(Phi*T^3)));
y=((1+mu+tau^2*rho_lin*gamma/(1/K*trace(Phi*T^2)))*1/rho_lin)/((1-tau^2)*(1+mu)+tau^2*mu*(csi+delta)^2/csi^2);
