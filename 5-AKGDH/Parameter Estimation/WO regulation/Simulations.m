clear all
close all
clc
%Driver
%load data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(nmol/min)
A_data = Data(:,2); %C_ICITm (M)
B_data = Data(:,3); %C_NADm (M)
C_data = Data(:,4); %C_NADHm (M) 

mpar = load('mpar.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
AKG_tot1=[0:0.001:0.6];
COA_tot1=[0.005 0.01 0.02 0.05];
NAD_tot1=[0.033 0.066 0.133 0.333];
for i=1:1:length(COA_tot1)
    for j=1:1:length(AKG_tot1)
    A_tot = AKG_tot1(j)/10^3; %C_AKGm (M)
    B_tot = COA_tot1(i)/10^3; %C_COAm (M)
    C_tot = NAD_tot1(i)/10^3; %C_NADm (M)
    Conc=[A_tot, B_tot, C_tot,0,0];
    flux1(i,j)=2*10^3*flux(mpar,Conc);
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 900 600])
subplot(2,3,1)
plot(A_data(1:4,:),Flux_Data(1:4,:),'k*',A_data(5:9,:),Flux_Data(5:9,:),'ko',A_data(10:14,:),Flux_Data(10:14,:),'k<',A_data(15:18,:),Flux_Data(15:18,:),'kd',AKG_tot1,flux1(1,:),'k',AKG_tot1,flux1(2,:),'k',AKG_tot1,flux1(3,:),'k',AKG_tot1,flux1(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('AKG (mM)')
ylabel('AKGDH Forward Flux(\mumol/min)')
legend('5 COA,33 NAD (\muM)','10 COA,66 NAD (\muM)','20 COA,133 NAD (\muM)','50 COA,333 NAD (\muM)','Location','north')
text(0.02,1.2,'A','Fontsize',18)
legend box off
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
COA_tot2=[0:0.001:0.1];
AKG_tot2=[0.025 0.05 0.1 0.5];
NAD_tot2=[0.02 0.04 0.08 0.4];
for i=1:1:length(AKG_tot2)
    for j=1:1:length(COA_tot2)
    A_tot = AKG_tot2(i)/10^3; %C_AKGm (M)
    B_tot = COA_tot2(j)/10^3; %C_COAm (M)
    C_tot = NAD_tot2(i)/10^3; %C_NADm (M)
    Conc=[A_tot, B_tot, C_tot,0,0];
    flux2(i,j)=10^3*flux(mpar,Conc);
    end
end
subplot(2,3,2)
plot(B_data(19:23,:),Flux_Data(19:23,:),'k*',B_data(24:28,:),Flux_Data(24:28,:),'ko',B_data(29:33,:),Flux_Data(29:33,:),'k<',B_data(34:38,:),Flux_Data(34:38,:),'kd',COA_tot2,flux2(1,:),'k',COA_tot2,flux2(2,:),'k',COA_tot2,flux2(3,:),'k',COA_tot2,1.6*flux2(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('COA (mM)')
legend('0.025 AKG,0.02 NAD (mM)','0.05 AKG,0.04 NAD (mM)','0.1 AKG,0.08 NAD (mM)','0.5 AKG,0.4 NAD (mM)','Location','north')
text(0.01,1,'B','Fontsize',18)
legend box off
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3 
NAD_tot3=[0:0.001:0.4];
AKG_tot3=[0.05 0.1 0.2 0.5];
COA_tot3=[0.005 0.01 0.02 0.05];
for i=1:1:length(AKG_tot3)
    for j=1:1:length(NAD_tot3)
    A_tot = AKG_tot3(i)/10^3; %C_AKGm (M)
    B_tot = COA_tot3(i)/10^3; %C_COAm (M)
    C_tot = NAD_tot3(j)/10^3; %C_NADm (M)
    Conc=[A_tot, B_tot, C_tot,0,0];
    flux3(i,j)=1.45*10^3*flux(mpar,Conc);
    end
end
subplot(2,3,3)
plot(C_data(39:43,:),Flux_Data(39:43,:),'k*',C_data(44:48,:),Flux_Data(44:48,:),'k<',C_data(49:53,:),Flux_Data(49:53,:),'ko',C_data(54:58,:),Flux_Data(54:58,:),'kd',NAD_tot3,flux3(1,:),'k',NAD_tot3,flux3(2,:),'k',NAD_tot3,flux3(3,:),'k',NAD_tot3,flux3(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NAD (mM)')
legend('0.05 AKG,0.005 COA (mM)','0.1 AKG,0.01 COA (mM)','0.2 AKG,0.02 COA (mM)','0.5 AKG,0.05 COA (mM)','Location','north')
text(0.02,1,'C','Fontsize',18)
legend box off
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 4 
NAD_tot4=[0:0.001:0.4];
NADH_tot4=[0 0.01 0.02 0.05];
for i=1:1:length(NADH_tot4)
    for j=1:1:length(NAD_tot4)
    A_tot = 0.5/10^3; %C_AKGm (M)
    B_tot = 0.5/10^3; %C_COAm (M)
    C_tot = NAD_tot4(j)/10^3; %C_NADm (M)
    E_tot = NADH_tot4(i)/10^3; %C_NADHm (M)
    
    Conc=[A_tot, B_tot, C_tot,0,E_tot];
    flux4(i,j)=10^3*flux(mpar,Conc);
    end
end
subplot(2,3,4)
plot(C_data(59:62,:),Flux_Data(59:62,:),'k*',C_data(63:67,:),Flux_Data(63:67,:),'ko',C_data(68:71,:),Flux_Data(68:71,:),'k<',C_data(72:75,:),Flux_Data(72:75,:),'kd',NAD_tot4,1.4*flux4(1,:),'k',NAD_tot4,1.2*flux4(2,:),'k',NAD_tot4,flux4(3,:),'k',NAD_tot4,0.8*flux4(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NAD (mM)')
ylabel('AKGDH Forward Flux (\mumol//min)')
legend('0 NADH(mM)','0.01 NADH(mM)','0.02 NADH(mM)','0.05 NADH(mM)','Location','southeast')
text(0.02,0.8,'D','Fontsize',18)
legend box off
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 5 
COA_tot5=[0:0.0001:0.01];
SCOA_tot5=[0 6.6 13.2];
for i=1:1:length(SCOA_tot5)
    for j=1:1:length(COA_tot5)
    A_tot = 1/10^3; %C_AKGm (M)
    B_tot = COA_tot5(j)/10^3; %C_COAm (M)
    C_tot = 0.336/10^3; %C_NADm (M)
    D_tot = SCOA_tot5(i)/10^6; %C_SCOAm (M)
    
    Conc=[A_tot, B_tot, C_tot,D_tot,0];
    flux5(i,j)=0.7*10^3*flux(mpar,Conc);
    end
end
subplot(2,3,5)
plot(B_data(76:79,:),Flux_Data(76:79,:),'k*',B_data(80:83,:),Flux_Data(80:83,:),'ko',B_data(84:87,:),Flux_Data(84:87,:),'k<',COA_tot5,flux5(1,:),'k',COA_tot5,flux5(2,:),'k',COA_tot5,flux5(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('COA (mM)')
legend('0 SCOA(\muM)','6.6 SCOA(\muM)','13.2 SCOA(\muM)','Location','southeast')
text(0.0005,0.3,'E','Fontsize',18)
legend box off
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 6 
COA_tot6=[0:0.0001:0.04];
SCOA_tot6=[0 6.6 13.2];
for i=1:1:length(SCOA_tot6)
    for j=1:1:length(COA_tot6)
    A_tot = 1/10^3; %C_AKGm (M)
    B_tot = COA_tot6(j)/10^3; %C_COAm (M)
    C_tot = 0.016/10^3; %C_NADm (M)
    D_tot = SCOA_tot6(i)/10^6; %C_NADm (M)
    E_tot = 0.042/10^3; %C_NADm (M)
    
    Conc=[A_tot, B_tot, C_tot,D_tot,E_tot];
    flux6(i,j)=3*10^3*flux(mpar,Conc);
    end
end
subplot(2,3,6)
plot(B_data(88:93,:),Flux_Data(88:93,:),'k*',B_data(94:99,:),Flux_Data(94:99,:),'k<',B_data(100:105,:),Flux_Data(100:105,:),'ko',COA_tot6,flux6(1,:),'k',COA_tot6,flux6(2,:),'k',COA_tot6,flux6(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('COA (mM)')
%ylabel('AKGDH Forward Flux (\mumol//min/mg protein)')
ylim([0 0.6])
legend('0 SCOA(\muM)','6.6 SCOA(\muM)','13.2 SCOA(\muM)','Location','southeast')
text(0.002,0.6,'F','Fontsize',18)
legend box off
box off
