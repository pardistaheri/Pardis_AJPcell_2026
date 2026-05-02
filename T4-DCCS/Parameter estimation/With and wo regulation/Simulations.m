%%%Load experimental data
Data = load('Data.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min)
B_data = Data(:,2); %C_Pie (M)

%%%Load Parameters
mpar = load('mpar.txt');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 
B_tot = [0:0.001:1]; %Suc
E_tot=[0 0.5];
for i=1:1:length(E_tot)
    for j=1:1:length(B_tot)
    B = B_tot(j)/10^3;
    E = E_tot(i)/10^3;
   
    Conc=[B, E];
    flux2(i,j)=flux(mpar,Conc)/10^3;
    end
end

h1 = figure(1);
set(h1,'Position',[10 10 300 300])
plot(B_data(7:11,:),Flux_Data(7:11,:),'k*',B_data(12:16,:),Flux_Data(12:16,:),'ko',...
    B_tot,flux2(1,:),'k',B_tot,flux2(2,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('SUC_e (mM)')
ylabel('Influx of SUC(mmol/mg)')
legend('Malonate_e = 0mM','Malonate_e = 0.5mM','Location','northwest')
%text(0.05,0.05,'B','Fontsize',18)
box off
legend box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure 1
% B_tot = [0:0.01:1]; %Cae
% for i=1:1:length(B_tot)
%     B = B_tot(i)/10^3;
%     E = 0;
%     Conc=[B, E];
%     flux1(:,i)=flux(mpar,Conc)/10^3;
% end
% 
% h1 = figure(1);
% set(h1,'Position',[10 10 600 300])
% subplot(1,2,1)
% plot(B_data(1:6,:),Flux_Data(1:6,:),'k*',B_tot,0.9*flux1,'k','LineWidth',1.5,'MarkerSize',8.0);
% set(gcf,'color','w')
% set(gca,'FontSize',14,'LineWidth',2.0)
% xlabel('SUC_e (mM)')
% ylabel('Influx of SUC(mmol/mg)')
% text(0.05,0.03,'A','Fontsize',18)
% box off

