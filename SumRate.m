function C = SumRate(H,a,sigma2,S,Pow)
updateSchedule = ["power","uniform","aa"];
numIterations = 100;
numRealizations = 20;
beta=1/4;%large-scale fading
[Ms,K,~,~] = size(H);
c = zeros(K,S);
b=2;
[V_RKA,us] = functionRKA(Ms,K,Pow,numRealizations,numIterations,H(:,:,1,1),updateSchedule(2));
[V_RKA2,us2] = functionRKA(Ms,K,Pow,numRealizations,numIterations,H(:,:,1,1),updateSchedule(1));
for s1=1:1:S
    for idx1 = 1:1:K
        if a==1
            P=RZF(H(:,:,1,1)',Pow);
        elseif a==2
            % [V_RKA,us] = functionRKA(Ms,K,Pow,numRealizations,numIterations(1,b,2),H(:,:,s1,s1)',updateSchedule(2));
            aa=reshape(V_RKA(:,3,:),[Ms K]);
            P=sqrt(Pow/trace(aa*aa'))*aa;
        else
            % [V_RKA2,us2] = functionRKA(Ms,K,Pow,numRealizations,numIterations(1,b,2),H(:,:,s1,s1)',updateSchedule(1));
            aa2=reshape(V_RKA2(:,3,:),[Ms K]);
            P=sqrt(Pow/trace(aa2*aa2'))*aa2;
        end
        ds = abs(H(:,idx1,1,1)'*P(:,idx1))^2;
        int = 0;
        for s2=1:1:S
            for idx2 = 1:1:K
                if a==1
                    P=RZF(H(:,:,1,1)',Pow);

                elseif a==2
                    %  [V_RKA,us] = functionRKA(Ms,K,Pow,numRealizations,numIterations(1,b,2),H(:,:,s1,s2)',updateSchedule(2));
                    aa=reshape(V_RKA(:,3,:),[Ms K]);
                    P=sqrt(Pow/trace(aa*aa'))*aa;
                else
                    %[V_RKA2,us2] = functionRKA(Ms,K,Pow,numRealizations,numIterations(1,b,2),H(:,:,s1,s2)',updateSchedule(1));
                    aa2=reshape(V_RKA2(:,3,:),[Ms K]);
                    P=sqrt(Pow/trace(aa2*aa2'))*aa2;
                end
                if s1==s2

                    int = int + abs(H(:,idx1,1,1)'*P(:,idx2))^2;
                else
                    int = int + beta*abs(H(:,idx1,1,1)'*P(:,idx2))^2;
                end
            end
        end
        sinr_k = ds/(sigma2+int-ds);
        c(idx1,s1) = log2(1+sinr_k);
    end
end
C = sum(sum(c))/(S);