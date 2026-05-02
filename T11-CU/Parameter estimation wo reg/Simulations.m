% CI biochemical reaction model fit and simulation
% the flux for CU reaction in this model was devolped in Audi & Dash lab previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2013
% experimental condition : temp = 298.15K, pH = 8, deltaSi=0, 
% The units are as follows Concentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg, 5 parameters
%%% Compute the free cation and metabolite concentrations and dye fluorescence
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1); %Flux(nM/min/ug)
A_data = Data(:,2); %C_Hm (M)

A_tot = [0:0.1:40]; %Cae

for i=1:1:length(A_tot)
    A = A_tot(i)/10^6;
    C = 0;
    E = 0;
    Conc=[A, C, E];
    flux1(i,:)=Model(mpar,Conc)*10^6;
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
%subplot(1,3,1);
plot(A_tot,flux1*1.05,'k',A_data(:,:),Flux_Data(:,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Ca_e (\muM)')
ylabel('Influx (nmol/mg)')
text(2,0.35,'A','Fontsize',14)
% title('Mg_e=Pi_e=0 mM')
box off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure 2 
% Data = load('Data_B.txt'); 
% Flux_Data = Data(:,1); %Flux(nM/min/ug)
% A_data = Data(:,2); %C_Hm (M)
% 
% A_tot = [0:0.1:40]; %Cae
% C_tot = [0 1];
% for i=1:1:length(C_tot)
%     for j=1:1:length(A_tot)
%         A = A_tot(j)/10^6;
%         C = C_tot(i)/10^3;
%         E = 1/10^3;
%         Conc=[A, C, E];
%         flux2(i,j)=Model(mpar,Conc)*10^6;
% end
% 
% subplot(1,3,2);
% plot(A_tot,flux2(1,:),'r',A_tot,flux2(2,:),'k',A_data(1:20,:),Flux_Data(1:20,:),'ro',A_data(21:41,:),Flux_Data(21:41,:),'k*','LineWidth',2.0,'MarkerSize',10.0);
% %plot(D,flux1(1,:),'r',D,flux1(2,:),'b','LineWidth',2.0,'MarkerSize',10.0);
% set(gcf,'color','w')
% set(gca,'FontSize',16)
% xlabel('Ca_e (\muM)')
% ylabel('Influx (nmol/mg)')
% title('Pi_e=1 mM')
% legend('Mg_e=0mM','Mg_e=1mM')
% legend box off
% box off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure 3 
% Data = load('Data_C.txt'); 
% Flux_Data = Data(:,1); %Flux(nM/min/ug)
% A_data = Data(:,2); %C_Hm (M)
% 
% A_tot = [0:0.1:200]; %Cae
% C_tot = [2 5];
% E_tot = [5 2];
% for i=1:1:length(C_tot)
%     for j=1:1:length(A_tot)
%         A = A_tot(j)/10^6;
%         C = C_tot(i)/10^3;
%         E = E_tot(i)/10^3;
%         Conc=[A, C, E];
%         flux3(i,j)=Model(mpar,Conc)*10^6;
% end
% 
% subplot(1,3,3);
% plot(A_tot,flux3(1,:),'r',A_tot,flux3(2,:),'k',A_data(1:17,:),Flux_Data(1:17,:),'ro',A_data(18:25,:),Flux_Data(18:25,:),'k*','LineWidth',2.0,'MarkerSize',10.0);
% %plot(D,flux1(1,:),'r',D,flux1(2,:),'b','LineWidth',2.0,'MarkerSize',10.0);
% set(gcf,'color','w')
% set(gca,'FontSize',16)
% xlabel('Ca_e (\muM)')
% ylabel('Influx (nmol/mg)')
% legend('Mg_e=2mM,Pi_e=5mM','Mg_e=1mM,Pi_e=2mM')
% legend box off
% box off