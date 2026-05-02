% Load experimental data
%Driver

mpar = load('mpar.txt');
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1); %Flux(muM/min)
B = Data(:,3); %C_NADm (M)

NAD_tot=[0:0.01:0.6]; %
MAL_tot=[0.33 0.5 0.67 1 2 4]; %

for i=1:1:length(MAL_tot)
    for j=1:1:length(NAD_tot)
        A_tot = MAL_tot(i)/10^3;
        B_tot = NAD_tot(j)/10^3;
        C_tot = 0;
        D_tot = 0;
    Conc=[A_tot, B_tot,C_tot, D_tot];
    flux1(i,j)=flux(mpar,Conc)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 900 600])
subplot(2,3,1)
plot(B(1:6,:),Flux_Data(1:6,:),'k*',B(7:12,:),Flux_Data(7:12,:),'ko',B(13:17,:),Flux_Data(13:17,:),'k<',B(18:23,:),Flux_Data(18:23,:),'kd',B(24:29,:),Flux_Data(24:29,:),'kp',B(30:35,:),Flux_Data(30:35,:),'k>',...
    NAD_tot,flux1(1,:),'k',NAD_tot,flux1(2,:),'k',NAD_tot,flux1(3,:),'k',NAD_tot,flux1(4,:),'k',NAD_tot,flux1(5,:),'k',NAD_tot,flux1(6,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NAD(mM)')
ylabel('MDH Forward Flux (\mumol/min/\mug)')
legend('0.33 mM MAL','0.5 mM MAL','0.67 mM MAL','1 mM MAL','2 mM MAL','4 mM MAL','Location','northwest')
ylim([0 0.2])
text(0.03,0.2,'A','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 
Data = load('Data_B.txt'); 
Flux_Data = Data(:,1); %Flux(muM/min)
A = Data(:,2); %C_MALm (M)

MAL_tot=[0:0.1:5]; %
NADH_tot=[0.005 0.00375 0.0025 0.00125 0]; %

for i=1:1:length(NADH_tot)
    for j=1:1:length(MAL_tot)
        A_tot = MAL_tot(j)/10^3;
        B_tot = 0.2/10^3;
        C_tot = 0;
        D_tot = NADH_tot(i)/10^3;
    Conc=[A_tot, B_tot,C_tot, D_tot];
    flux2(i,j)=0.9*flux(mpar,Conc)*10^3;
    end
end

subplot(2,3,2)
plot(A(1:6,:),Flux_Data(1:6,:),'k*',A(7:12,:),Flux_Data(7:12,:),'ko',A(13:18,:),Flux_Data(13:18,:),'k<',A(19:24,:),Flux_Data(19:24,:),'kd',A(25:30,:),Flux_Data(25:30,:),'kp',...
    MAL_tot,flux2(1,:),'k',MAL_tot,flux2(2,:),'k',MAL_tot,flux2(3,:),'k',MAL_tot,flux2(4,:),'k',MAL_tot,flux2(5,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('MAL(mM)')
%ylabel('Forward Flux (\mumol/min/\mug)')
legend('5 \muM NADH','3.75 \muM NADH','2.5 \muM NADH','1.25 \muM NADH','0 \muM NADH','Location','northwest')
ylim([0 0.1])
text(0.3,0.1,'B','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 3 
Data = load('Data_C.txt'); 
Flux_Data = Data(:,1); %Flux(muM/min)
A = Data(:,2); %C_MALm (M)

MAL_tot=[0:0.1:5]; %
OXA_tot=[0.008 0.006 0.004 0.002 0]; %

for i=1:1:length(OXA_tot)
    for j=1:1:length(MAL_tot)
        A_tot = MAL_tot(j)/10^3;
        B_tot = 0.25/10^3;
        C_tot = OXA_tot(i)/10^3;
        D_tot = 0;
    Conc=[A_tot, B_tot,C_tot, D_tot];
    flux3(i,j)=flux(mpar,Conc)*10^3;
    end
end
subplot(2,3,3)
plot(A(1:6,:),Flux_Data(1:6,:),'k*',A(7:12,:),Flux_Data(7:12,:),'ko',A(13:18,:),Flux_Data(13:18,:),'k<',A(19:24,:),Flux_Data(19:24,:),'kd',A(25:30,:),Flux_Data(25:30,:),'kp',...
    MAL_tot,flux3(1,:),'k',MAL_tot,flux3(2,:),'k',MAL_tot,flux3(3,:),'k',MAL_tot,flux3(4,:),'k',MAL_tot,flux3(5,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('MAL(mM)')
%ylabel('Forward Flux (\mumol/min/\mug)')
legend('8 \muM OXA','6 \muM OXA','4 \muM OXA','2 \muM OXA','0 \muM OXA','Location','northwest')
ylim([0 0.2])
text(0.3,0.2,'C','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 4 
Data = load('Data_D.txt'); 
Flux_Data = Data(:,1); %Flux(muM/min)
B = Data(:,3); %C_NADm (M)

NAD_tot=[0:0.01:0.6]; %
NADH_tot=[0.00375 0.0025 0.00125 0]; %

for i=1:1:length(NADH_tot)
    for j=1:1:length(NAD_tot)
        A_tot = 2/10^3;
        B_tot = NAD_tot(j)/10^3;
        C_tot = 0;
        D_tot = NADH_tot(i)/10^3;
    Conc=[A_tot, B_tot,C_tot, D_tot];
    flux4(i,j)=0.55*flux(mpar,Conc)*10^3;
    end
end
subplot(2,3,4)
plot(B(1:6,:),Flux_Data(1:6,:),'k*',B(7:12,:),Flux_Data(7:12,:),'ko',B(13:18,:),Flux_Data(13:18,:),'k<',B(19:24,:),Flux_Data(19:24,:),'kd',...
    NAD_tot,flux4(1,:),'k',NAD_tot,flux4(2,:),'k',NAD_tot,flux4(3,:),'k',NAD_tot,flux4(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NAD(mM)')
ylabel('MDH Forward Flux (\mumol/min/\mug)')
legend('3.75 \muM NADH','2.5 \muM NADH','1.25 \muM NADH','0 \muM NADH','Location','northwest')
ylim([0 0.1])
text(0.03,0.1,'D','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 5 
Data = load('Data_E.txt'); 
Flux_Data = Data(:,1); %Flux(muM/min)
B = Data(:,3); %C_NADm (M)

NAD_tot=[0:0.01:0.6]; %
OXA_tot=[0.006 0.0045 0.003 0.0015 0]; %

for i=1:1:length(OXA_tot)
    for j=1:1:length(NAD_tot)
        A_tot = 1/10^3;
        B_tot = NAD_tot(j)/10^3;
        C_tot = OXA_tot(i)/10^3;
        D_tot = 0;
    Conc=[A_tot, B_tot,C_tot, D_tot];
    flux5(i,j)=0.9*flux(mpar,Conc)*10^3;
    end
end

subplot(2,3,5)
plot(B(1:6,:),Flux_Data(1:6,:),'k*',B(7:12,:),Flux_Data(7:12,:),'ko',B(13:18,:),Flux_Data(13:18,:),'k<',B(19:24,:),Flux_Data(19:24,:),'kd',B(25:30,:),Flux_Data(25:30,:),'kp',...
    NAD_tot,flux5(1,:),'k',NAD_tot,flux5(2,:),'k',NAD_tot,flux5(3,:),'k',NAD_tot,flux5(4,:),'k',NAD_tot,1.1*flux5(5,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NAD(mM)')
%ylabel('Forward Flux (\mumol/min/\mug)')
legend('6 \muM OXA','4.5 \muM OXA','3 \muM OXA','1.5 \muM OXA','0 \muM OXA','Location','northwest')
ylim([0 0.15])
text(0.03,0.15,'E','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 6 
Data = load('Data_Rev.txt'); 
Flux_rev = Data(:,1); %Flux(muM/min)
C = Data(:,4); %C_OXAm (M)

OXA_tot=[0:0.001:0.06]; %
NADH_tot=[0.0067 0.01 0.015 0.02 0.025 0.05]; %

for i=1:1:length(NADH_tot)
    for j=1:1:length(OXA_tot)
        A_tot = 0;
        B_tot = 0;
        C_tot = OXA_tot(j)/10^3;
        D_tot = NADH_tot(i)/10^3;
    Conc=[A_tot, B_tot,C_tot, D_tot];
    flux6(i,j)=-flux(mpar,Conc)*10^15;
    end
end
subplot(2,3,6)
plot(C(1:7,:),Flux_rev(1:7,:),'k*',C(8:14,:),Flux_rev(8:14,:),'ko',C(15:21,:),Flux_rev(15:21,:),'k<',C(22:28,:),Flux_rev(22:28,:),'kd',C(29:35,:),Flux_rev(29:35,:),'kp',C(36:42,:),Flux_rev(36:42,:),'k>',OXA_tot,1.7*flux6(1,:),'k',OXA_tot,1.7*flux6(2,:),'k',OXA_tot,2.3*flux6(3,:),'k',OXA_tot,3*flux6(4,:),'k',OXA_tot,3.5*flux6(5,:),'k',OXA_tot,5.2*flux6(6,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('OXA(mM)')
ylabel('MDH Reverse Flux (\mumol/min/\mug)')
legend('0.0067 mM NADH','0.01 mM NADH','0.015 mM NADH','0.02 mM NADH','0.025 mM NADH','0.05 mM NADH','Location','northwest')
ylim([0 2])
text(0.005,2,'F','Fontsize',18)
legend boxoff
box off