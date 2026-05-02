clear all
close all
clc
%Driver
Data = load('Data.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min)
A_Data = Data(:,2); %C_ACOAm (\muM)
B_Data = Data(:,3); %C_OXAm (\muM)
C_Data = Data(:,4); %C_ACOAm (\muM)
D_Data = Data(:,5); %C_CITm (mM)

mpar = load('mpar.txt');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 CIT=COA=0
OXA_tot=[2.5 5 10];
ACOA_tot=[0:0.1:50];
for i=1:1:length(OXA_tot)
    for j=1:1:length(ACOA_tot)
        A_tot = ACOA_tot(j)/10^6;
        B_tot = OXA_tot(i)/10^6;
        C_tot = 0;
        D_tot = 0;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux1(i,j)=flux(mpar,Conc);
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 1200 300])
subplot(1,4,1);
plot(A_Data(1:8,:),Flux_Data(1:8,:),'k*',A_Data(9:16,:),Flux_Data(9:16,:),'ko',A_Data(17:23,:),Flux_Data(17:23,:),'k<',ACOA_tot,flux1(1,:),'k',ACOA_tot,flux1(2,:),'k',ACOA_tot,flux1(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('ACOA (\muM)')
ylabel('CITS Forward Flux (mmol/min)')
legend('2.5 \muM OXA','5 \muM OXA','10 \muM OXA','Location','southeast')
text(4,0.16,'A','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 CIT=COA=0
ACOA_tot=[2.5 5 10];
OXA_tot=[0:0.1:35];
for i=1:1:length(ACOA_tot)
    for j=1:1:length(OXA_tot)
        A_tot = ACOA_tot(i)/10^6;
        B_tot = OXA_tot(j)/10^6;
        C_tot = 0;
        D_tot = 0;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux2(i,j)=flux(mpar,Conc);
    end
end
subplot(1,4,2);
plot(B_Data(24:31,:),Flux_Data(24:31,:),'k*',B_Data(32:39,:),Flux_Data(32:39,:),'ko',B_Data(40:47,:),Flux_Data(40:47,:),'k<',OXA_tot,flux2(1,:),'k',OXA_tot,flux2(2,:),'k',OXA_tot,flux2(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('OXA (\muM)')
%ylabel('CITS Forward Flux (mmol/min)')
legend('2.5 \muM ACOA','5 \muM ACOA','10 \muM ACOA','Location','southeast')
text(2,0.15,'B','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3 ACOA=OXA=0
CIT_tot=[1 2 4 10];
COA_tot=[0:1:100];
for i=1:1:length(CIT_tot)
    for j=1:1:length(COA_tot)
        A_tot = 0;
        B_tot = 0;
        C_tot = COA_tot(j)/10^6;
        D_tot = CIT_tot(i)/10^3;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux3(i,j)=-flux(mpar,Conc)*10^3;
    end
end

subplot(1,4,3);
plot(C_Data(48:51,:),-Flux_Data(48:51,:)*10^3,'k*',C_Data(52:55,:),-Flux_Data(52:55,:)*10^3,'ko',C_Data(56:59,:),-Flux_Data(56:59,:)*10^3,'k<',C_Data(60:63,:),-Flux_Data(60:63,:)*10^3,'kd',COA_tot,7*flux3(1,:),'k',COA_tot,7*flux3(2,:),'k',COA_tot,7*flux3(3,:),'k',COA_tot,8*flux3(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('COA (\muM)')
ylabel('CITS Reverse Flux (\mumol/min)')
legend('1 mM CIT','2 mM CIT','4 mM CIT','10 mM CIT','Location','northwest')
text(5,0.05,'C','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 4 ACoA=OXA=0
COA_tot=[20 33 50 100];
CIT_tot=[0:0.01:11];
for i=1:1:length(COA_tot)
    for j=1:1:length(CIT_tot)
        A_tot = 0;
        B_tot = 0;
        C_tot = COA_tot(i)/10^6;
        D_tot = CIT_tot(j)/10^3;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux4(i,j)=-flux(mpar,Conc)*10^3;
    end
end
subplot(1,4,4);
plot(D_Data(64:67,:),-Flux_Data(64:67,:)*10^3,'k*',D_Data(68:71,:),-Flux_Data(68:71,:)*10^3,'ko',D_Data(72:75,:),-Flux_Data(72:75,:)*10^3,'k<',D_Data(76:79,:),-Flux_Data(76:79,:)*10^3,'kd',CIT_tot,7*flux4(1,:),'k',CIT_tot,7*flux4(2,:),'k',CIT_tot,7*flux4(3,:),'k',CIT_tot,7*flux4(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('CIT (mM)')
%ylabel('CITS Reverse Flux (mmol/min)')
legend('20 \muM COA','33 \muM COA','50 \muM COA','100 \muM COA','Location','southeast')
text(0.8,0.05,'D','Fontsize',18)
legend boxoff
box off
