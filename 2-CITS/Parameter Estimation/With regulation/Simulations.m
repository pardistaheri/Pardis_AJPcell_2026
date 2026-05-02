clear all
close all
clc


Data = load('Data2.txt'); 
Flux_Data = Data(:,1)/10^1; %Flux(mM/min)
A_Data = Data(:,2); %C_ACOAm (\muM)
B_Data = Data(:,3); %C_OXAm (\muM)
ATP_Data = Data(:,4); %C_ATPm (mM)
ADP_Data = Data(:,5); %C_ADPm (mM)
AMP_Data = Data(:,6); %C_AMPm (mM)
SCOA_Data = Data(:,7); %C_SCOAm (\muM)

mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
ACOA_tot=[0:0.1:60];
OXA_tot = [1.3 2.5 3.5 5.6 10.6];
for i=1:1:length(OXA_tot)
    for j=1:1:length(ACOA_tot)
        A = ACOA_tot(j)/10^6;
        B = OXA_tot(i)/10^6;
        ATP = 0; %C_ATPm (M)
        ADP = 0; %C_ADPm (M)
        AMP = 0; %C_AMPm (M)
        SCOA = 0; %C_SCOAm (M)
    Conc=[A, B, ATP,ADP,AMP,SCOA];
    flux1(i,j)=flux(mpar,Conc);
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 900 600])
subplot(2,3,1);
plot(A_Data(1:5,:),Flux_Data(1:5,:),'k*',A_Data(6:10,:),Flux_Data(6:10,:),'ko',A_Data(11:15,:),Flux_Data(11:15,:),'k<',A_Data(16:20,:),Flux_Data(16:20,:),'kd',A_Data(21:25,:),Flux_Data(21:25,:),'kp',ACOA_tot,flux1(1,:),'k',ACOA_tot,flux1(2,:),'k',ACOA_tot,flux1(3,:),'k',ACOA_tot,flux1(4,:),'k',ACOA_tot,flux1(5,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ACOA (\muM)')
ylabel('CITS Forward Flux (arbitrary unit)')
legend('1.3 mM OXA','2.5 mM OXA','3.5 mM OXA','5.6 mM OXA','10.6 mM OXA','Location','northwest')
text(4,0.5,'E','Fontsize',18)
ylim([0 0.5])
legend boxoff
box off
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
ACOA_tot=[0:0.1:80];
ATP_tot = [0 4 0 0];
ADP_tot = [0 0 4 0];
AMP_tot = [0 0 0 4];
for i=1:1:length(ATP_tot)
    for j=1:1:length(ACOA_tot)
        A = ACOA_tot(j)/10^6;
        B = 34/10^6;
        ATP = ATP_tot(i)/10^3; %C_ATPm (M)
        ADP = ADP_tot(i)/10^3; %C_ADPm (M)
        AMP = AMP_tot(i)/10^3; %C_AMPm (M)
        SCOA = 0; %C_SCOAm (M)
    Conc=[A, B, ATP,ADP,AMP,SCOA];
    flux2(i,j)=flux(mpar,Conc);
    end
end
subplot(2,3,2);
plot(A_Data(26:30,:),Flux_Data(26:30,:),'k*',A_Data(31:35,:),Flux_Data(31:35,:),'ko',A_Data(36:40,:),Flux_Data(36:40,:),'k<',A_Data(41:45,:),Flux_Data(41:45,:),'kd',ACOA_tot,0.9*flux2(1,:),'k',ACOA_tot,flux2(2,:),'k',ACOA_tot,flux2(3,:),'k',ACOA_tot,flux2(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ACOA (\muM)')
%ylabel('CITS Forward Flux (arbitrary unit)')
legend('0 mM AXP','4 mM ATP','4 mM ADP','4 mM AMP','Location','northwest')
text(6,0.5,'F','Fontsize',18)
ylim([0 0.5])
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3 
OXA_tot=[0:0.1:22];
ATP_tot = [0 4 0 0];
ADP_tot = [0 0 4 0];
AMP_tot = [0 0 0 4];
for i=1:1:length(ATP_tot)
    for j=1:1:length(OXA_tot)
        A = 84/10^6;
        B = OXA_tot(j)/10^6;
        ATP = ATP_tot(i)/10^3; %C_ATPm (M)
        ADP = ADP_tot(i)/10^3; %C_ADPm (M)
        AMP = AMP_tot(i)/10^3; %C_AMPm (M)
        SCOA = 0; %C_SCOAm (M)
    Conc=[A, B, ATP,ADP,AMP,SCOA];
    flux3(i,j)=flux(mpar,Conc);
    end
end
subplot(2,3,3);
plot(B_Data(46:50,:),Flux_Data(46:50,:),'k*',B_Data(51:55,:),Flux_Data(51:55,:),'ko',B_Data(56:60,:),Flux_Data(56:60,:),'k<',B_Data(61:65,:),Flux_Data(61:65,:),'kd',OXA_tot,0.75*flux3(1,:),'k',OXA_tot,1.2*flux3(2,:),'k',OXA_tot,1.2*flux3(3,:),'k',OXA_tot,flux3(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('OXA (\muM)')
%ylabel('CITS Forward Flux (arbitrary unit)')
legend('0 mM AXP','4 mM ATP','4 mM ADP','4 mM AMP','Location','northwest')
text(1.5,0.5,'G','Fontsize',18)
ylim([0 0.5])
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 4 
ATP_tot=[0:0.01:4.5];
OXA_tot = [0.7 2.8 10.3];
for i=1:1:length(OXA_tot)
    for j=1:1:length(ATP_tot)
        A = 60/10^6;
        B = OXA_tot(i)/10^6;
        ATP = ATP_tot(j)/10^3; %C_ATPm (M)
        ADP = 0; %C_ADPm (M)
        AMP = 0; %C_AMPm (M)
        SCOA = 0; %C_SCOAm (M)
    Conc=[A, B, ATP,ADP,AMP,SCOA];
    flux4(i,j)=flux(mpar,Conc);
    end
end
subplot(2,3,4);
plot(ATP_Data(66:70,:),Flux_Data(66:70,:),'k*',ATP_Data(71:75,:),Flux_Data(71:75,:),'ko',ATP_Data(76:80,:),Flux_Data(76:80,:),'k<',ATP_tot,flux4(1,:),'k',ATP_tot,1.5*flux4(2,:),'k',ATP_tot,1.8*flux4(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ATP (mM)')
ylabel('CITS Forward Flux (arbitrary unit)')
legend('0.7 \muM OXA','2.8 \muM OXA','10.3 \muM OXA','Location','northeast')
text(0.22,0.5,'H','Fontsize',18)
ylim([0 0.5])
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 5 
ACOA_tot=[0:0.1:40];
SCOA_tot = [0 80 160];
for i=1:1:length(SCOA_tot)
    for j=1:1:length(ACOA_tot)
        A = ACOA_tot(j)/10^6;
        B = 3.1/10^6;
        ATP = 0; %C_ATPm (M)
        ADP = 0; %C_ADPm (M)
        AMP = 0; %C_AMPm (M)
        SCOA = SCOA_tot(i)/10^6; %C_SCOAm (M)
    Conc=[A, B, ATP,ADP,AMP,SCOA];
    flux5(i,j)=flux(mpar,Conc);
    end
end
subplot(2,3,5);
plot(A_Data(81:86,:),Flux_Data(81:86,:),'k*',A_Data(87:92,:),Flux_Data(87:92,:),'ko',A_Data(93:98,:),Flux_Data(93:98,:),'k<',ACOA_tot,flux5(1,:),'k',ACOA_tot,flux5(2,:),'k',ACOA_tot,0.95*flux5(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ACOA (\muM)')
%ylabel('CITS Forward Flux (arbitrary unit)')
legend('0 \muM SCOA','80 \muM SCOA','160 \muM SCOA','Location','southeast')
text(4,0.2,'I','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 6 
mpar_pH = load('mpar_pH.txt');
%figure 1 pH regulation
Data_pH = load('Data_pH.txt'); 
Flux_Data = Data_pH(:,1); %Flux(arbitrary unit)
pH_data = Data_pH(:,2); %pHm (M)

pH_tot=[5:0.01:9];
for i=1:1:length(pH_tot)
    pH = pH_tot(i);
    flux1(i,:)=flux_pH(mpar_pH,pH);
end
subplot(2,3,6);
plot(pH_tot,flux1,'k',pH_data(:,:),Flux_Data(:,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('CITS Forward Flux (normalized)')
text(5.2,1,'J','Fontsize',18)
box off

