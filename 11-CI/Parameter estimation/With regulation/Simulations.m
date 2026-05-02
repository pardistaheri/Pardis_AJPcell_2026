% CI biochemical reaction model fit and simulation for Ca2+ regulations
% Model fitted to data from Sadek et al. 2004
% experimental condition : temp = 298.15K, deltaSi=0, 
% The units are as follows Concentration M, Flux mmol/nin, Volume =1 ml, Mass in Microg, 2 parameters

mpar = load('mpar.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1 Ca regulation
Data = load('Data_Ca.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min/ug)
Ca_data = Data(:,3); %C_Cam (M)

Ca=[0:0.01:100];
for i=1:1:length(Ca)
    Ca_tot = Ca(i)/10^6;
    A_tot = 20/10^6;
    Conc=[Ca_tot,A_tot];
    flux2(i,:)=flux(mpar,Conc)*10^3;
end
h1 = figure(1)
set(h1,'Position',[10 10 900 300])
subplot(1,3,1);
plot(Ca,flux2,'k',Ca_data(1:7,:),Flux_Data(1:7,:),'kd','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Ca^2^+ (\muM)')
ylabel('CI Forward Flux (\mumol/min)')
text(5,1.1,'F','Fontsize',18)
ylim([0 1.1])
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 NADH & Ca regulation
Data = load('Data_Ca.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min/ug)
A_data = Data(:,2); %C_Cam (M)

Ca=[0 7];
A = [0:0.1:20];
for i=1:1:length(Ca)
    for j=1:1:length(A)
    Ca_tot = Ca(i)/10^6;
    A_tot = A(j)/10^6;
    Conc=[Ca_tot,A_tot];
    flux2(i,j)=flux(mpar,Conc)*10^3;
    end
end
subplot(1,3,2);
plot(A,flux2(1,:),'k',A,flux2(2,:),'--k',A_data(8:13,:),Flux_Data(8:13,:),'ko',A_data(14:19,:),Flux_Data(14:19,:),'kd','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
legend('0 \muM Ca^2^+','20 \muM Ca^2^+','Location','southeast')
xlabel('NADH (\muM)')
%ylabel('Flux (\mumol/min)')
text(1,1.1,'G','Fontsize',18)
ylim([0 1.1])
legend box off
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure 3 Ca simulation
% Ca=[0 4 8 12 16 20]/10^6;
% A = [0:0.1:20]/10^6;
% for i=1:1:length(Ca)
%     for j=1:1:length(A)
%     Ca_tot = Ca(i);
%     A_tot = A(j);
%     Conc=[Ca_tot,A_tot];
%     flux2(i,j)=flux(mpar,Conc)*10^3;
%     end
% subplot(1,3,3);    
% plot(A*10^6,flux2(i,:),'LineWidth',2.0); hold on;
% set(gcf,'color','w')
% set(gca,'FontSize',12)
% legend('0 \muM Ca^2^+','4 \muM Ca^2^+','8 \muM Ca^2^+','12 \muM Ca^2^+','16 \muM Ca^2^+','20 \muM Ca^2^+','Location','southeast')
% xlabel('NADH (\muM)')
% ylabel('Flux (\mumol/min)')
% ylim([0 1.1])
% legend box off
% box off
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpar = load('mpar_pH.txt');
%figure 1 pH regulation
Data = load('Data_pH.txt'); 
Flux_Data = Data(:,1)/10^3; %Flux(mM/min/ug)
pH_data = Data(:,2); %pHm (M)
pH=[6:0.01:9];
for i=1:1:length(pH)
    pH_tot = pH(i);
    Conc=[pH_tot, 0];
    flux4(i,:)=flux_pH(mpar,Conc)*10^3;
end

subplot(1,3,3);   
plot(pH,flux4,'k',pH_data(1:5,:),Flux_Data(1:5,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
%ylabel('Flux (\mumol/min)')
text(6.2,1.2,'H','Fontsize',18)
ylim([0.2 1.2])
box off