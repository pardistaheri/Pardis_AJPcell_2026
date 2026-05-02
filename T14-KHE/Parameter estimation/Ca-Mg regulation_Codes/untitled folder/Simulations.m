% 1. CI biochemical reaction model fit
% the flux for CITS reaction in this model was devolped in Audi & Dash lab
% previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2010
% The units are as follows Cncentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg,5 parameters


%Driver
mpar = load('mpar.txt');

Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux%
Mg_data = Data(:,3); %Mg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1 
Mg=[10:0.01:1000];
for i=1:1:length(Mg)
    Mg_tot = Mg(i)/10^6;
    flux2(i,:)=Model(mpar,Mg_tot);
end
plot(log10(Mg),flux2,'k',log10(Mg_data(7:13,:)),Flux_Data(7:13,:),'kd','LineWidth',2.5,'MarkerSize',10.0);
set(gcf,'color','w')
set(gca,'FontSize',12)
xlabel('logMg(\muM)')
ylabel('%Activity')
box off

