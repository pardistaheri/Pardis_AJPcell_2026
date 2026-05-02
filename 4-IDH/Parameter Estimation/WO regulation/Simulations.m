clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(nmol/min)
A_data = Data(:,2); %C_ICITm (M)
B_data = Data(:,3); %C_NADm (M)
D_data = Data(:,4); %C_NADHm (M) 

mpar = load('mpar.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
ICIT_tot=[0:0.01:2.5];
for i=1:1:length(ICIT_tot)
    A_tot = ICIT_tot(i)/10^3; %C_ICITm (M)
    B_tot = 2/10^3; %C_NADm (M)
    D_tot = 0; %C_NADHm (M)
    Conc=[A_tot, B_tot, D_tot,0,0];
    flux1(:,i)=10^7*flux(mpar,Conc);
end
h1 = figure(1)
set(h1,'Position',[10 10 900 300])
subplot(1,3,1)
plot(ICIT_tot,0.9*flux1,'k',A_data(1:14,:),Flux_Data(1:14,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ICIT (mM)')
ylabel('IDH Forward Flux (munits/mg protein)')
text(0.1,60,'A','Fontsize',18)
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
NADH_tot = [0 39];
NAD_tot=[0:0.01:0.8];
for i=1:1:length(NADH_tot)
    for j=1:1:length (NAD_tot)
    A_tot = 0.45/1000; %C_ICITm (M)
    B_tot = NAD_tot(j)/10^3; %C_NADm (M)
    D_tot = NADH_tot(i)/10^6; %C_NADHm (M)
    Conc=[A_tot, B_tot,D_tot,0,0];
    flux2(i,j)=10^6*flux(mpar,Conc);
    end
end
subplot(1,3,2)
plot(B_data(15:25,:),Flux_Data(15:25,:),'k*',B_data(26:33,:),Flux_Data(26:33,:),'ko',NAD_tot,flux2(1,:),'k',NAD_tot,flux2(2,:)*2,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NAD (mM)')
%ylabel('IDH Forward Flux (munits/mg protein)')
legend('NADH=0 mM','NADH= 0.039 mM','Location','northwest')
text(0.02,1.2,'B','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3
NADH_tot=[0:0.1:50];
for i=1:1:length(NADH_tot)
    A_tot = 0.45/1000; %C_ICITm (M)
    B_tot = 0.25/1000; %C_NADm (M)
    D_tot = NADH_tot(i)/10^6; %C_NADHm (M)
    Conc=[A_tot, B_tot,D_tot,0,0];
    flux3(:,i)=10^6*flux(mpar,Conc);
end
subplot(1,3,3)
plot(NADH_tot,flux3,'k',D_data(34:43,:),Flux_Data(34:43,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NADH (\muM)')
%ylabel('IDH Forward Flux (munits/mg protein)')
text(2,0.4,'C','Fontsize',18)
box off


