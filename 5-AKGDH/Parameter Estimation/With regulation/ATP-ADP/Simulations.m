clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(nmol/min)
A_data = Data(:,2); %C_ICITm (M)
ADP_data = Data(:,5); %C_NADm (M)
ATP_data = Data(:,6); %C_NADHm (M) 

mpar = load('mpar.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
ATP_tot = [0 1.5];
ADP_tot=[0:0.01:5];
for i=1:1:length(ATP_tot)
    for j=1:1:length (ADP_tot)
    A_tot = 2/10^3; %C_AKGDHm (M)
    B_tot = 0.25/10^3; %C_COAm (M)
    C_tot = 0.1/10^3; %C_NADm (M)
    ADP = ADP_tot(j)/10^3;
    ATP = ATP_tot(i)/10^3;
    Conc=[A_tot, B_tot,C_tot,ADP,ATP];
    flux1(i,j)=flux(mpar,Conc);
    end
end

ATP_tot2=[0:0.01:5];
for j=1:1:length (ATP_tot2)
    A_tot = 2/10^3; %C_AKGDHm (M)
    B_tot = 0.25/10^3; %C_COAm (M)
    C_tot = 0.1/10^3; %C_NADm (M)
    ADP = 0/10^3;
    ATP = ATP_tot2(j)/10^3;
    Conc=[A_tot, B_tot,C_tot,ADP,ATP];
    flux2(:,j)=flux(mpar,Conc);
end    
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1);
plot(ADP_data(1:8,:),Flux_Data(1:8,:),'k*',ADP_data(9:16,:),Flux_Data(9:16,:),'ko',ATP_data(50:56,:),Flux_Data(50:56,:),'k<',ADP_tot,0.8*flux1(1,:),'k',ADP_tot,1.1*flux1(2,:),'k',ATP_tot2,flux2(:,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ADP or ATP (mM)')
ylabel('AKGDH Forward Flux (unit/mg protein)')
legend('ADP with 0mM ATP','ADP with 1.5 mM ATP','ATP','Location','northeast')
ylim([0 120])
text(0.5,120,'G','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
ATP_tot3 = [0 0.6 0];
ADP_tot3 = [0 0 0.6];
AKG_tot3 = [0:0.01:11];

for i=1:1:length(ATP_tot3)
    for j=1:1:length (AKG_tot3)
    A_tot = AKG_tot3(j)/10^3; %C_AKGDHm (M)
    B_tot = 0.25/10^3; %C_COAm (M)
    C_tot = 1/10^3; %C_NADm (M)
    ADP = ADP_tot3(i)/10^3;
    ATP = ATP_tot3(i)/10^3;
    Conc=[A_tot, B_tot,C_tot,ADP,ATP];
    flux3(i,j)=1.2*flux(mpar,Conc);
    end
end
    
subplot(1,2,2);
plot(A_data(17:27,:),Flux_Data(17:27,:),'k*',A_data(28:38,:),Flux_Data(28:38,:),'ko',A_data(39:49,:),Flux_Data(39:49,:),'k<',AKG_tot3,flux3(1,:),'k',AKG_tot3,flux3(2,:),'k',AKG_tot3,flux3(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('AKG (mM)')
ylabel('AKGDH Forward Flux (nmol/mg protein)')
legend('ADP=ATP=0 mM','ATP=0.6 mM','ADP=0.6 mM','Location','southeast')
text(1,140,'H','Fontsize',18)
legend boxoff
box off
