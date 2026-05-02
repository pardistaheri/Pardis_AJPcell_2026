% CI biochemical reaction model fit and simulation
% the flux for CI reaction in this model was devolped in Audi & Dash lab previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2013
% experimental condition : temp = 298.15K, pH = 8, deltaSi=0, 
% The units are as follows Concentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg, 5 parameters
 
%%% Compute the free cation and metabolite concentrations and dye fluorescence
mpar = load('mpar.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1 
Data = load('Data_C.txt'); 
Flux_Data = Data(:,1); %Flux(1/s)
A = Data(:,2); %C_UQ (M)

dPsi_tot = [126 175 178 166];
pHi_tot = [7.35 6.5 7.35 8.35];
fract_tot=[0:0.001:0.8];

for i=1:1:length(dPsi_tot)
    for j=1:1:length(fract_tot)
        fract = fract_tot(j);
        B_tot = 200/10^6;
        dPsi = dPsi_tot(i);
        dpH = 0.4;
        pHi = pHi_tot(i);
        Conc=[fract, B_tot];
        flux3(i,j)=flux(mpar,dPsi,Conc,pHi,dpH)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
plot(A(1:13,:),Flux_Data(1:13,:),'k*',A(14:26,:),Flux_Data(14:26,:),'ko',A(27:40,:),Flux_Data(27:40,:),'k<',A(41:53,:),Flux_Data(41:53,:),'kd',fract_tot(1:180),flux3(1,1:180)*4,'k',fract_tot(1:500),flux3(2,1:500),'k',fract_tot,flux3(3,:),'k',fract_tot,flux3(4,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Fractional Cytochrome reduced')
ylabel('CIV Forward Flux (1/sec)')
legend('deltaPsi=126mV, pHi=7.35','deltaPsi=175mV, pHi=6.5','deltaPsi=178mV, pHi=7.35','deltaPsi=166mV, pHi=8.35','Location','northeast')
ylim([0 60]);
xlim([0 1]);
text(0.05,60,'C','Fontsize',18)
legend boxoff
box off