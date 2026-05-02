% KHE biochemical transporter model fit and simulation
% the flux for KHE transporter in this model was devolped in Audi & Dash lab previously
% Model fitted to data from Bazil et al. 2010
% experimental condition : temp = 315.15K
% The units are as follows Concentration M, Flux mmol/nin, Volume =1 ml, Mass in Microg
clear all
close all
clc


mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux(nmol/min/ug)
D_data = Data(:,4); %C_Cae (M)

C = [15 30]; %Ke
D = [0:0.01:100]; %Hm

for i=1:1:length(C)
    for j=1:1:length(D)
        B_tot = 10^-7.4;
        C_tot = C(i)/10^3;
        D_tot = D(j)/10^9;
        
        Conc=[0,B_tot, C_tot, D_tot];
        flux1(i,j)=Model(mpar,Conc)*10^6;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(D_data(1:8,:),Flux_Data(1:8,:),'k*',D_data(9:15,:),Flux_Data(9:15,:),'ko',D,flux1(1,:),'k',D,flux1(2,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
legend('K_e = 15 mM', 'K_e = 30 mM');
xlabel('H_m(nM)')
ylabel('Flux (nmol/min/mg)')
text(5,40,'A','Fontsize',16)
legend box off
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure 2
% Data = load('Data.txt'); 
% Flux_Data = Data(:,1); %Flux(nmol/min/ug)
% D_data = Data(:,4); %C_Cae (M)
% 
% D = [30:1:178]; %Hm
% B = [9.05:-(9.05-7.6)/length(D):7.6]; %He
% 
% for i=1:1:length(D)
%     A_tot = 0;
%     B_tot = 10^-B(i);
%     C_tot = 20/10^3;
%     D_tot = D(i)/10^9;
%         
%         Conc=[0,B_tot, C_tot, D_tot];
%         flux2(:,i)=Model(mpar,Conc)*10^6;
% end
% subplot(1,2,2)
% plot(D,flux2(1,:),'k',D_data(16:20,:),Flux_Data(16:20,:),'k*','LineWidth',2.0,'MarkerSize',10.0);
% %plot(D,flux2(1,:),'k','LineWidth',2.0,'MarkerSize',10.0);
% set(gcf,'color','w')
% set(gca,'FontSize',16)
% xlabel('H_m(nM)')
% ylabel('Flux (nmol/min/mg)')
% box off