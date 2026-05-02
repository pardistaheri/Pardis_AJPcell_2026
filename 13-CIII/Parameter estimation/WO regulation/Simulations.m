%%% Compute the free cation and metabolite concentrations and dye fluorescence
mpar = load('mpar.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1 CytCr=UQ=0
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min/ug)
A = Data(:,2); %C_CytCo (M)

CytCo_tot=[0:0.01:48];
UQH2_tot=[25 15 10 7]/10^6;
CytCr_tot=0;
for i=1:1:length(UQH2_tot)
    for j=1:1:length(CytCo_tot)
        A_tot = CytCo_tot(j)/10^6;
        B_tot = UQH2_tot(i);
        C_tot = CytCr_tot;
        D_tot = 0;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux1(i,j)=flux(mpar,Conc)/10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 900 300])
subplot(1,3,1)
plot(A(1:13,:),Flux_Data(1:13,:),'k*',A(14:27,:),Flux_Data(14:27,:),'ko',A(28:40,:),Flux_Data(28:40,:),'k<',A(41:53,:),Flux_Data(41:53,:),'kd',...
    CytCo_tot,flux1(1,:),'k',CytCo_tot,flux1(2,:),'k',CytCo_tot,flux1(3,:),'k',CytCo_tot,flux1(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('CytCo (\muM)')
ylabel('CIII Forward Flux (mmol/sec)')
legend('25 \muM UQH2','15 \muM UQH2','10 \muM UQH2','7 \muM UQH2','Location','southeast')
ylim([0 0.5])
xlim([0 48])
text(1,0.5,'A','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 
Data = load('Data_B.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux
B = Data(:,3); %C_UQ (M)

CytCo_tot = [30 20 12 7]/10^6;
UQH2_tot=[0:0.01:32];
CytCr_tot=10/10^6;
for i=1:1:length(CytCo_tot)
    for j=1:1:length(UQH2_tot)
        A_tot = CytCo_tot(i);
        B_tot = UQH2_tot(j)/10^6;
        C_tot = CytCr_tot;
        D_tot = 20*10^-6;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux4(i,j)=0.9*flux(mpar,Conc)/10^3;
    end
end
subplot(1,3,2);
plot(B(1:12,:),Flux_Data(1:12,:),'k*',B(13:24,:),Flux_Data(13:24,:),'ko',B(25:37,:),Flux_Data(25:37,:),'k<',B(38:50,:),Flux_Data(38:50,:),'kd',...
    UQH2_tot,flux4(1,:),'k',UQH2_tot,flux4(2,:),'k',UQH2_tot,flux4(3,:),'k',UQH2_tot,flux4(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('UQH_2 (\muM)')
% ylabel('CIII Forward Flux (mmol/sec)')
legend('30 \muM CytCo','20 \muM CytCo','12 \muM CytCo','7 \muM CytCo','Location','southeast')
text(1,0.25,'B','Fontsize',18)
ylim([0 0.25])
xlim([0 32])
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 3 
Data = load('Data_C.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min/ug)
A = Data(:,2); %C_CytCo (M)

CytCo_tot=[0:0.01:38];
UQH2_tot=20/10^6;
CytCr_tot=[0 10]/10^6;
for i=1:1:length(CytCr_tot)
    for j=1:1:length(CytCo_tot)
        A_tot = CytCo_tot(j)/10^6;
        B_tot = UQH2_tot;
        C_tot = CytCr_tot(i);
        D_tot = 0*10^-6;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux2(i,j)=1.15*flux(mpar,Conc)/10^3;
    end
end
subplot(1,3,3);
plot(A(1:9,:),Flux_Data(1:9,:),'k*',A(10:20,:),Flux_Data(10:20,:),'k<',...
    CytCo_tot,flux2(1,:),'k',CytCo_tot,flux2(2,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('CytCo (\muM)')
% ylabel('CIII Forward Flux (mmol/sec)')
legend('0 \muM CytCr','10 \muM CytCr','Location','southeast')
text(1,0.5,'C','Fontsize',18)
ylim([0 0.5])
xlim([0 32])
legend boxoff
box off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 4 
Data = load('Data_D.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min/ug)
B = Data(:,3); %C_UQH2 (M)

UQH2_tot=[0:0.01:30];
CytCo_tot=10/10^6;
CytCr_tot=[0 5 10 15]/10^6;
for i=1:1:length(CytCr_tot)
    for j=1:1:length(UQH2_tot)
        A_tot = CytCo_tot;
        B_tot = UQH2_tot(j)/10^6;
        C_tot = CytCr_tot(i);
        D_tot = 5/10^6;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux3(i,j)=1.15*flux(mpar,Conc)/10^3;
    end
end
h2 = figure(2)
set(h2,'Position',[10 10 300 300])
plot(B(1:12,:),Flux_Data(1:12,:),'k*',B(13:24,:),Flux_Data(13:24,:),'ko',B(25:36,:),Flux_Data(25:36,:),'k<',B(37:48,:),Flux_Data(37:48,:),'kd',...
    UQH2_tot,flux3(1,:),'k',UQH2_tot,flux3(2,:),'k',UQH2_tot,flux3(3,:),'k',UQH2_tot,flux3(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('UQH_2 (\muM)')
ylabel('CIII Forward Flux (mmol/sec)')
legend('0 \muM CytCr','5 \muM CytCr','10 \muM CytCr','15 \muM CytCr','Location','southeast')
text(1,0.4,'D','Fontsize',18)
ylim([0 0.4])
xlim([0 30])
legend boxoff
box off
