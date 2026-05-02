%%%Load experimental data
Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux(1/s)
A_Data = Data(:,2)/10^3; %C_ADPe (M)
D_Data = Data(:,3)/10^3; %C_ATPe (M)

%%%Load Parameters
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1 
A_tot=[0:0.001:0.11]; %ADPe uM
D_tot = [0 100 400]; %ATPe uM

for i=1:1:length(D_tot)
    for j=1:1:length(A_tot)
        A = A_tot(j)/10^3;
        B = 5/10^3;
        C = 5/10^3;
        D = D_tot(i)/10^6;
        dPsi = 180;
        Conc=[A,B,C,D,dPsi];
        flux180(i,j)=flux(mpar,Conc);
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 900 300])
subplot(1,3,1)
plot(A_Data(1:6,:),Flux_Data(1:6,:),'k*',A_Data(7:12,:),Flux_Data(7:12,:),'ko',A_Data(13:18,:),Flux_Data(13:18,:),'k<',...
    A_tot,flux180(1,:),'k',A_tot,flux180(2,:),'k',A_tot,flux180(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
%title('dPsi=180mV');
xlabel('ADP_e mM')
ylabel('Influx of ADP(1/S)')
legend('0 \muM ATP_e','100 \muM ATP_e', '400 \muM ATP_e','Location','southeast');
ylim([0 30])
text(0.005,30,'A','Fontsize',18)
legend boxoff
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2
A_tot=[0:0.001:0.11]; %ADP_o uM
D_tot =[0 20 100]; %ATP_o uM

for i=1:1:length(D_tot)
    for j=1:1:length(A_tot)
        A = A_tot(j)/10^3;
        B = 5/10^3;
        C = 5/10^3;
        D = D_tot(i)/10^6;
        dPsi = 0;
        Conc=[A,B,C,D,dPsi];
        flux0(i,j)=flux(mpar,Conc);
    end
end

subplot(1,3,2)
plot(A_Data(19:24,:),Flux_Data(19:24,:),'k*',A_Data(25:30,:),Flux_Data(25:30,:),'ko',A_Data(31:36,:),Flux_Data(31:36,:),'k<',...
    A_tot,flux0(1,:),'k',A_tot,flux0(2,:),'k',A_tot,flux0(3,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
%title('dPsi=0mV');
xlabel('ADP_e mM')
%ylabel('Flux (1/S)')
legend('0 \muM ATP_e','20 \muM ATP_e', '100 \muM ATP_e','Location','southeast');
ylim([0 12])
text(0.005,12,'B','Fontsize',18);
legend boxoff
box off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 3 

A_tot = [0:0.01:0.3]; %ADP_o uM
D_tot = [0:0.01:0.3]; %ATP_o uM
dPsi_tot = [180 0];
for i=1:1:length(dPsi_tot)
    for j=1:1:length(A_tot)
        A = A_tot(j)/10^3;
        B = 5/10^3;
        C = 5/10^3;
        D = D_tot(j)/10^3;
        dPsi = dPsi_tot(i);
        Conc=[A,B,C,D,dPsi];
        flux1(i,j)=1.06*flux(mpar,Conc);
    end
end
subplot(1,3,3)
plot(A_Data(37:41,:),Flux_Data(37:41,:),'k*',A_Data(42:45,:),Flux_Data(42:45,:),'ko',...
    A_tot,flux1(1,:),'k',A_tot,flux1(2,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
%title('ADP=ATP')
xlabel('ADP_e mM')
%ylabel('Flux (1/S)')
legend('180 mV dPsi','0 mV dPsi','Location','northwest');
ylim([0 35])
text(0.01,35,'C','Fontsize',18);
legend boxoff
box off