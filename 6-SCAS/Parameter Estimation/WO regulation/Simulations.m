clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(nmol/min)
A_data = Data(:,2); %C_SCOAm (M)
B_data = Data(:,3); %C_GDPm (M)
E_data = Data(:,6); %C_GTPm (M) 

mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
GDP_tot=[0:0.001:0.1];
GTP_tot=[0 0.05 0.15];
for i=1:1:length(GTP_tot)
    for j=1:1:length(GDP_tot)
    A_tot = 0.1/10^3; %C_SCOA (M)
    B_tot = GDP_tot(j)/10^3; %C_GDP (M)
    C_tot = 50/10^3; %C_Pi (M)
    D_tot = 0; %C_SUCm (M)
    E_tot = GTP_tot(i)/10^3; %C_GTPm (M)
    F_tot = 0; %C_COAm (M)
    Conc=[A_tot, B_tot, C_tot,D_tot, E_tot, F_tot];
    flux1(i,j)=flux(mpar,Conc);
    end
end

h1 = figure(1)
set(h1,'Position',[10 10 900 300])
subplot(1,3,1)
plot(B_data(24:28,:),Flux_Data(24:28,:),'k*',B_data(29:32,:),Flux_Data(29:32,:),'ko',B_data(33:37,:),Flux_Data(33:37,:),'k<',...
    GDP_tot,flux1(1,:),'k',GDP_tot,flux1(2,:),'k',GDP_tot,flux1(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('GDP (mM)')
ylabel('SCAS Forward Flux (unitless)')
legend('0 mM GTP','0.05 mM GTP','0.15 mM GTP','Location','southeast')
text(0.01,0.04,'A','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
SCOA_tot=[0:0.001:0.18];
COA_tot=[0 0.02 0.04];
for i=1:1:length(COA_tot)
    for j=1:1:length(SCOA_tot)
    A_tot = SCOA_tot(j)/10^3; %C_SCOA (M)
    B_tot = 0.05/10^3; %C_GDP (M)
    C_tot = 1/10^3; %C_Pi (M)
    D_tot = 0; %C_SUCm (M)
    E_tot = 0; %C_GTPm (M)
    F_tot = COA_tot(i)/10^3; %C_COAm (M)
    Conc=[A_tot, B_tot, C_tot,D_tot, E_tot, F_tot];
    flux2(i,j)=flux(mpar,Conc);
    end
end
subplot(1,3,2)
plot(A_data(38:45,:),Flux_Data(38:45,:),'k*',A_data(46:51,:),Flux_Data(46:51,:),'ko',A_data(52:57,:),Flux_Data(52:57,:),'k<',...
    SCOA_tot,flux2(1,:),'k',SCOA_tot,flux2(2,:),'k',SCOA_tot,flux2(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('SCOA (mM)')
%ylabel('SCAS Forward Flux (unitless)')
legend('0 mM COA','0.02 mM COA','0.04 mM COA','Location','southeast')
text(0.01,0.015,'B','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3 
SCOA_tot=[0:0.001:0.18];
SUC_tot=[0 0.5 2];
for i=1:1:length(SUC_tot)
    for j=1:1:length(SCOA_tot)
    A_tot = SCOA_tot(j)/10^3; %C_SCOA (M)
    B_tot = 0.05/10^3; %C_GDP (M)
    C_tot = 1/10^3; %C_Pi (M)
    D_tot = SUC_tot(i)/10^3; %C_SUCm (M)
    E_tot = 0; %C_GTPm (M)
    F_tot = 0; %C_COAm (M)
    Conc=[A_tot, B_tot, C_tot,D_tot, E_tot, F_tot];
    flux3(i,j)=flux(mpar,Conc);
    end
end
subplot(1,3,3)
plot(A_data(58:65,:),Flux_Data(58:65,:),'k*',A_data(66:72,:),Flux_Data(66:72,:),'ko',A_data(73:79,:),Flux_Data(73:79,:),'k<',...
    SCOA_tot,flux3(1,:),'k',SCOA_tot,flux3(2,:),'k',SCOA_tot,flux3(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('SCOA (mM)')
%ylabel('SCAS Forward Flux (unitless)')
legend('0 mM SUC','0.5 mM SUC','2 mM SUC','Location','southeast')
text(0.01,0.0151,'C','Fontsize',18)
legend boxoff
box off