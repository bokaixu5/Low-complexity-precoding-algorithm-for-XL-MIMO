function delta=deterministic(Phi,t,K)
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

Dphi=eig(Phi);
delta=1/K*sum(Dphi);
M=length(Dphi);
delta_prev=0;
while(abs(delta-delta_prev)>1e-8)
delta_prev=delta;
delta=1/K*sum(Dphi./(1+t*Dphi/(1+t*delta_prev)));
end


