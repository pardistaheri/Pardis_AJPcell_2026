% Load experimental data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(mM/min)
A = Data(:,2); %C_FUMm (uM)

mpar = load('mpar.txt');
options = [];       % ODE options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
FUM_tot=[0:0.001:0.05]; %
ATP=[0 42 83 166]*10^-6; %

for i=1:1:length(ATP)
    for j=1:1:length(FUM_tot)
        A_tot = FUM_tot(j)*10^-3;
        ATP_tot = ATP(i);
    Conc=[A_tot, ATP_tot];
    flux1(i,j)=Flux(mpar,Conc)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1)
plot(A(1:3,:),Flux_Data(1:3,:),'k*',A(4:6,:),Flux_Data(4:6,:),'ko',A(7:9,:),Flux_Data(7:9,:),'k<',A(10:12,:),Flux_Data(10:12,:),'kd',...
    FUM_tot,flux1(1,:),'k',FUM_tot,flux1(2,:),'k',FUM_tot,flux1(3,:),'k',FUM_tot,flux1(4,:),'k', 'LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('FUM (mM)')
ylabel('FH Forward Flux (\mumol/min)')
legend('0 mM ATP', '0.042 mM ATP','0.083 mM ATP','0.166 mM ATP','Location','southeast')
text(0.004,4,'C','Fontsize',18)
legend boxoff
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure 2
mpar_pH = load('mpar_pH.txt');

Data_pH = load('Data_pH.txt'); 
Flux_Data = Data_pH(:,1); %Flux(arbitrary unit)
pH_data = Data_pH(:,2); %pHm (M)

pH_tot=[4.5:0.01:8.5];
for i=1:1:length(pH_tot)
    pH = pH_tot(i);
    flux1(i,:)=3.2*Flux_pH(mpar_pH,pH);
end
subplot(1,2,2)
plot(pH_tot,flux1,'k',pH_data(:,:),Flux_Data(:,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('pH')
ylabel('FH Forward Flux (normalized)')
text(4.2,1.2,'D','Fontsize',18)
box off