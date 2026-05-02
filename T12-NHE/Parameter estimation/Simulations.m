% CI biochemical reaction model fit and simulation
% the flux for CI reaction in this model was devolped in Audi & Dash lab previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2013
% experimental condition : temp = 298.15K, pH = 8, deltaSi=0, 
% The units are as follows Concentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg, 5 parameters
 
%%% Compute the free cation and metabolite concentrations and dye fluorescence
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
Data = load('Data_A.txt'); 
Flux_Data = Data(:,1); %Flux(nM/min/ug)
D_data = Data(:,2); %C_Hm (M)

C = [15 9.6]; %Nae
B = [50:0.1:500]; %He

for i=1:1:length(C)
    for j=1:1:length(B)
        B_tot = B(j)/10^9;
        C_tot = C(i)/10^3;
        D_tot = 100/10^9;
      
        Conc=[B_tot, C_tot, D_tot];
        flux1(i,j)=Model(mpar,Conc)*10^6;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1);
plot(D_data(1:8,:),Flux_Data(1:8,:),'k*',D_data(9:14,:),Flux_Data(9:14,:),'ko',B,flux1(1,:)*0.95,'k',B,flux1(2,:)*0.8,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
legend('Na_e = 15mM', 'Na_e = 9.6')
xlabel('H_e (nM)')
ylabel('Eflux (nmol/min)')
text(20,2.5,'A','Fontsize',16)
legend box off
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2 
Data = load('Data_B.txt'); 
Flux_Data = Data(:,1); 
C_data = Data(:,3); 

C=[0:0.01:50];
for i=1:1:length(C)
    B_tot = 100/10^9;
    C_tot = C(i)/10^3;
    D_tot = 112.2/10^9;
    Conc=[B_tot,C_tot,D_tot];
    flux2(:,i)=Model(mpar,Conc)*10^6;
end
subplot(1,2,2);
plot(C,flux2,'k',C_data(1:14,:),Flux_Data(1:14,:),'k*','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Na_e (mM)')
ylabel('Eflux (nmol/min)')
legend('H_e = 100 nM')
text(2,5,'B','Fontsize',16)
legend box off
box off
