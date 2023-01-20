load('ber3.mat','par','res');
marker_style = {'kx-','bo:','rs--','mv-.','gp-.','bs--','y*--'};
h = figure(1);
for d=1:length(par.precoder)
    semilogy(par.NTPdB_list,res.BER(d,:),marker_style{d},'LineWidth',2);
    if (d==1)
        hold on
    end
end
%hold off
grid on
box on
xlabel('Normalized transmit power [dB]','FontSize',12,'Interpreter','Latex')
ylabel('Uncoded bit error rate (BER)','FontSize',12,'Interpreter','Latex');
if length(par.NTPdB_list) > 1
    axis([min(par.NTPdB_list) max(par.NTPdB_list) 1e-3 1]);
end
legend(par.precoder,'FontSize',12,'Interpreter','Latex','location','northeast')
set(gca,'FontSize',12);
%if par.save
    % save eps figure (in color and with a reasonable bounding box)
%    print(h,'-loose','-depsc',[ par.simName '_' num2str(par.runId) ])
%end
clear;
load('ber1.mat','par','res');
marker_style = {'kx-','bo:','rs--','mv-.','gp-.','bs--','y*--'};
h = figure(1);
for d=1:length(par.precoder)
    semilogy(par.NTPdB_list,res.BER(d,:),marker_style{d},'LineWidth',2);
    if (d==1)
        hold on
    end
end
hold off
grid on
box on
xlabel('Normalized transmit power [dB]','FontSize',12,'Interpreter','Latex')
ylabel('Uncoded bit error rate (BER)','FontSize',12,'Interpreter','Latex');
if length(par.NTPdB_list) > 1
    axis([min(par.NTPdB_list) max(par.NTPdB_list) 1e-3 1]);
end
legend(par.precoder,'FontSize',12,'Interpreter','Latex','location','northeast')
set(gca,'FontSize',12);
%if par.save
    % save eps figure (in color and with a reasonable bounding box)
    %print(h,'-loose','-depsc',[ par.simName '_' num2str(par.runId) ])
%end
 set(gca,'xLim',[-10 40]);