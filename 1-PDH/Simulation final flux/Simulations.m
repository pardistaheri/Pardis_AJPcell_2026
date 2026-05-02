clear all
close all
clc
%Driver
Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux(uM/min/ug)
A = Data(:,2); %C_PYRm (mM)
B = Data(:,3); %C_COAm (mM)
C = Data(:,4); %C_NADm (mM)
D = 0; %C_ACOAm
E = 0; %C_CO2
F = 0; %C_NADHm
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 NAD constant
COA_tot=[0.2 0.1 0.07 0.05];
PYR_tot=[0:0.001:0.6];
for i=1:1:length(COA_tot)
    for j=1:1:length(PYR_tot)
        A_tot = PYR_tot(j)/1000;
        B_tot = COA_tot(i)/1000;
        C_tot = 0.05/1000;
    Conc=[A_tot, B_tot, C_tot];
    Pflux(i,j)=flux(Conc)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 900 300]);
subplot(1,3,1);
plot(A(1:4,:),Flux_Data(1:4,:),'k*',A(5:8,:),Flux_Data(5:8,:),'ko',A(9:12,:),Flux_Data(9:12,:),'k<',A(13:16,:),Flux_Data(13:16,:),'kd',...
    PYR_tot,Pflux(1,:),'k',PYR_tot,Pflux(2,:),'k',PYR_tot,Pflux(3,:),'k',PYR_tot,Pflux(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('PYR (mM)')
ylabel('PDH Forward Flux (\mumol/min)')
legend('0.20 mM COA','0.10 mM COA','0.07 mM COA','0.05 mM COA','Location','southeast')
text((0.6)*0.1,(0.05),'A','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 PYR constant
NAD_tot=[0.5 0.25 0.1 0.05];
COA_tot=[0:0.001:0.3];
PYR_tot=0.25;
for i=1:1:length(NAD_tot)
    for j=1:1:length(COA_tot)
        A_tot = PYR_tot/1000;
        B_tot = COA_tot(j)/1000;
        C_tot = NAD_tot(i)/1000;
    Conc=[A_tot, B_tot, C_tot];
    Cflux(i,j)=flux(Conc)*10^3;
    end
    
end
subplot(1,3,2);
plot(B(17:20,:),Flux_Data(17:20,:),'k*',B(21:24,:),Flux_Data(21:24,:),'ko',B(25:28,:),Flux_Data(25:28,:),'k<',B(29:32,:),Flux_Data(29:32,:),'kd',...
    COA_tot,Cflux(1,:),'k',COA_tot,Cflux(2,:),'k',COA_tot,Cflux(3,:),'k',COA_tot,Cflux(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('COA (mM)')
text((0.6)*0.05,(0.08),'B','Fontsize',18)
legend('0.5 mM NAD','0.25 mM NAD','0.1 mM NAD','0.05 mM NAD','Location','southeast')
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3 COA constant
PYR_tot=[0.5 0.25 0.1 0.05];
NAD_tot=[0:0.001:0.6];
COA_tot=0.1;
for i=1:1:length(PYR_tot)
    for j=1:1:length(NAD_tot)
        A_tot = PYR_tot(i)/1000;
        B_tot = COA_tot/1000;
        C_tot = NAD_tot(j)/1000;
    Conc=[A_tot, B_tot, C_tot];
    Nflux(i,j)=flux(Conc)*10^3;
    end
end
subplot(1,3,3);
plot(C(33:36,:),Flux_Data(33:36,:),'k*',C(37:40,:),Flux_Data(37:40,:),'ko',C(41:44,:),Flux_Data(41:44,:),'k<',C(45:48,:),Flux_Data(45:48,:),'kd',...
    NAD_tot,Nflux(1,:),'k',NAD_tot,Nflux(2,:),'k',NAD_tot,Nflux(3,:),'k',NAD_tot,Nflux(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NAD (mM)')
legend('0.5 mM PYR','0.25 mM PYR','0.1 mM PYR','0.05 mM PYR','Location','southeast')
text((0.6)*0.075,(0.08),'C','Fontsize',18)
legend boxoff
box off

