% Load experimental data
%Driver
Data_F = load('Data_F.txt'); 
Flux_Data = Data_F(:,1); %Flux(mM/min)
A_Data = Data_F(:,2); %C_FUMm (M)
B_data = Data_F(:,3); %C_MALm (M)

mpar_F = load('mpar_F.txt');
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 1 
FUM_tot=[0:0.001:0.4]; %

for i=1:1:length(FUM_tot)
    A_tot = FUM_tot(i)/10^3;
    B_tot = 0;
    Conc=[A_tot, B_tot];
    flux1(i,:)=Flux_F(mpar_F,Conc)*10^3;
    end

h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1)
plot(FUM_tot,flux1,'k',A_Data(1:4,:),Flux_Data(1:4,:),'ko','LineWidth',1.5,'MarkerSize',8);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('FUM (mM)')
ylabel('FH Forward Flux (\mumol/min)')
text(0.03,50,'A','Fontsize',18)
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data_R = load('Data_R.txt'); 
Flux_Data = Data_R(:,1); %Flux(mM/min)
A_Data = Data_R(:,2); %C_FUMm (M)
B_data = Data_R(:,3); %C_MALm (M)

mpar_R = load('mpar_R.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2
MAL_tot=[0:0.001:5]; %

for i=1:1:length(MAL_tot)
    A_tot = 0;
    B_tot = MAL_tot(i)/10^3;
    Conc=[A_tot, B_tot];
    flux2(i,:)=Flux_R(mpar_R,Conc)*10^3;
end
subplot(1,2,2)
plot(MAL_tot,flux2,'k',B_data(1:4,:),Flux_Data(1:4,:),'ko','LineWidth',1.5,'MarkerSize',8);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('MAL (mM)')
ylabel('FH Reverse Flux (\mumol/min)')
text(0.4,20,'B','Fontsize',18)
box off