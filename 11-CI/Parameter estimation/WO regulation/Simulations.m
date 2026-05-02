mpar = load('mpar.txt');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1 NAD=UQH2=0
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min/ug)
A = Data(:,2); %C_NADH (M)

NADH_tot=[0:0.01:12];
UQ_tot=[25 50 100 200]/10^6;
NAD_tot=0;
UQH2_tot=0;
for i=1:1:length(UQ_tot)
    for j=1:1:length(NADH_tot)
        A_tot = NADH_tot(j)/10^6;
        B_tot = UQ_tot(i);
        C_tot = NAD_tot;
        D_tot = UQH2_tot;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux1(i,j)=flux(mpar,Conc)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 900 600])
subplot(2,3,1);
plot(A(1:11,:),Flux_Data(1:11,:),'k*',A(12:21,:),Flux_Data(12:21,:),'ko',A(22:31,:),Flux_Data(22:31,:),'k<',A(32:39,:),Flux_Data(32:39,:),'kd',...
    NADH_tot,flux1(1,:),'k',NADH_tot,flux1(2,:),'k',NADH_tot,1.1*flux1(3,:),'k',NADH_tot,0.95*flux1(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NADH (\muM)')
ylabel('CI Forward Flux (\mumol/min)')
legend('25 \muM UQ','50 \muM UQ','100 \muM UQ','200 \muM UQ','Location','southeast')
text(1,1.3,'A','Fontsize',18)
ylim([0 1.3])
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 UQH2=0
Data = load('Data_B.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min/ug)
A = Data(:,2); %C_NADH (M)

NADH_tot=[0:0.01:12];
UQ_tot=100/10^6;
NAD_tot=[0 200 400 600]/10^6;
UQH2_tot=0;
for i=1:1:length(NAD_tot)
    for j=1:1:length(NADH_tot)
        A_tot = NADH_tot(j)/10^6;
        B_tot = UQ_tot;
        C_tot = NAD_tot(i);
        D_tot = UQH2_tot;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux2(i,j)=1.1*flux(mpar,Conc)*1000;
    end
end
subplot(2,3,2);
plot(A(1:10,:),Flux_Data(1:10,:),'k*',A(11:20,:),Flux_Data(11:20,:),'ko',A(21:30,:),Flux_Data(21:30,:),'k<',A(31:39,:),Flux_Data(31:39,:),'kd',...
    NADH_tot,flux2(1,:),'k',NADH_tot,flux2(2,:),'k',NADH_tot,flux2(3,:),'k',NADH_tot,flux2(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NADH (\muM)')
% ylabel('Flux \mumol/min')
legend('0 \muM NAD','200 \muM NAD','400 \muM NAD','600 \muM NAD','Location','southeast')
text(1,1.2,'B','Fontsize',18)
legend boxoff
box off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 3 NAD=0
Data = load('Data_C.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min/ug)
A = Data(:,2); %C_NADH (M)

NADH_tot=[0:0.01:12];
UQ_tot=100/10^6;
NAD_tot =0;
UQH2_tot=[0 100 200]/10^6;
for i=1:1:length(UQH2_tot)
    for j=1:1:length(NADH_tot)
        A_tot = NADH_tot(j)/10^6;
        B_tot = UQ_tot;
        C_tot = NAD_tot;
        D_tot = UQH2_tot(i);
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux3(i,j)=1.05*flux(mpar,Conc)*1000;
    end
end
subplot(2,3,3);
plot(A(1:10,:),Flux_Data(1:10,:),'k*',A(11:20,:),Flux_Data(11:20,:),'ko',A(21:30,:),Flux_Data(21:30,:),'k<',...
    NADH_tot,1.05*flux3(1,:),'k',NADH_tot,flux3(2,:),'k',NADH_tot,flux3(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('NADH (\muM)')
% ylabel('Flux \mumol/min')
legend('0 \muM UQH2','100 \muM UQH2','200 \muM UQH2','Location','southeast')
text(1,1.2,'C','Fontsize',18)
legend boxoff
box off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 4 UQH2=0
Data = load('Data_D.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min/ug)
B = Data(:,3); %C_UQ (M)

NADH_tot = 10/10^6;
UQ_tot=[0:0.1:220];
NAD_tot=[0 400 600]/10^6;
UQH2_tot =0;
for i=1:1:length(NAD_tot)
    for j=1:1:length(UQ_tot)
        A_tot = NADH_tot;
        B_tot = UQ_tot(j)/10^6;
        C_tot = NAD_tot(i);
        D_tot = UQH2_tot;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux4(i,j)=flux(mpar,Conc)*1000;
    end
end
subplot(2,3,4);
plot(B(1:6,:),Flux_Data(1:6,:),'k*',B(7:11,:),Flux_Data(7:11,:),'ko',B(12:17,:),Flux_Data(12:17,:),'k<',...
    UQ_tot,flux4(1,:),'k',UQ_tot,flux4(2,:),'k',UQ_tot,flux4(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('UQ (\muM)')
ylabel('CI Forward Flux \mumol/min')
legend('0 \muM NAD','400 \muM NAD','600 \muM NAD','Location','southeast')
text(10,1.4,'D','Fontsize',18)
ylim([0 1.4])
legend boxoff
box off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 5 NAD=0
Data = load('Data_E.txt'); 
Flux_Data = Data(:,1); %Flux(mM/min/ug)
B = Data(:,3); %C_UQ (M)

NADH_tot = 10/10^6;
UQH2_tot=[0 100 200]/10^6;
NAD_tot =0;
UQ_tot=[0:0.1:220];
for i=1:1:length(UQH2_tot)
    for j=1:1:length(UQ_tot)
        A_tot = NADH_tot;
        B_tot = UQ_tot(j)/10^6;
        C_tot = NAD_tot;
        D_tot = UQH2_tot(i);
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux5(i,j)=flux(mpar,Conc)*1000;
    end
end
subplot(2,3,5);
plot(B(1:6,:),Flux_Data(1:6,:),'k*',B(7:12,:),Flux_Data(7:12,:),'ko',B(13:17,:),Flux_Data(13:17,:),'k<',...
    UQ_tot,flux5(1,:),'k',UQ_tot,flux5(2,:),'k',UQ_tot,flux5(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('UQ (\muM)')
% ylabel('Flux \mumol/min')
legend('0 \muM UQH_2','100 \muM UQH_2','200 \muM UQH_2','Location','southeast')
text(10,1.4,'E','Fontsize',18)
ylim([0 1.4])
legend boxoff
box off
