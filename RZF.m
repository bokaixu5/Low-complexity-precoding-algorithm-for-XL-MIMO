function P = RZF(H,pow)
[K,M] = size(H);
pre = H'*inv(H*H'+ K/pow*eye(K));
P = sqrt(pow/trace(pre*pre'))*pre;