clear all
M = 256;
K = 10;
SNR = 10;
xi = 1 / SNR;
T = 200;
for i = 1 : 10
    Q = randn(M, K) + 1j * randn(M, K);
    MMSED = (Q' * Q + xi * eye(K))^(-1)  * Q';
    for j = 1 : M
        y = zeros(M,1);
        y(j) = 1;
        for r = 1 : K
            %p(r)=1/K;
            p(r) = (norm(Q(:,r))^2 + xi) / (norm(Q, 'fro') + K * xi);
        end
        b = Q' * y;
        QH = Q';
        p = p / sum(p);
        idx_set = 1 : K;
        u = zeros(M,1);
        v = zeros(K,1);
        for t = 1 : T
            r = randsrc(1,1, [idx_set;p]);
            qrH = QH(r,:);
            gamma = (b(r) - qrH * u -  xi * v(r)) / (norm(qrH)^2 + xi);
            u = u + gamma * qrH';
            v(r) = v(r) + gamma;
            KAD(:, j, t) = v;
        end
    end
    for t = 1 : T
        MSE(i, t) = norm(squeeze(KAD(:, :, t)) - MMSED, 'fro')^2;
    end
end
load('x_1','a')
figure
semilogy(1:T, mean(MSE),'R','LineWidth',2)
hold on
semilogy(1:T, a,'B-.','LineWidth',2)
set(gca,'xLim',[80 180]);
set(gca,'FontSize',12);
xlabel('Number of iterations','Interpreter','Latex')
ylabel('$\left \|x^{t} -x_{0}   \right \|^{2}_{2}  $','Interpreter','Latex')
legend('SwoR-rKA','rKA','Location','best','Interpreter','Latex')
grid on