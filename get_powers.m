%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Anum Ali 
%
% Last Modified: Feb, 2019
%
% If you use this code or any (modified) part of it in any publication,
% please cite the paper: Anum Ali, Elisabeth de Carvalho, and Robert W.
% Heath Jr., "Linear Receivers in Non-stationary Massive MIMO Channels with
% Visibility Regions", IEEE Wireless Communications Letters.
%
% Contact person email: anumali@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description:
% A function to get user powers to satisfy the constraint in equation (2) in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% K: The number of users
% P: The total transmit power (scalar)
% etas: diag(G'G) from equation (2) in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Argument:
% Pmtx: The KxK diagonal matrix with user powers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pmtx=get_powers(K,P,etas)
cvx_begin quiet
variables pvecc(K) epsilon
minimize epsilon
subject to
(etas'*pvecc-P)^2<=epsilon;
pvecc>=0;
cvx_end
Pmtx=diag(pvecc);
end