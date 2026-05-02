% 1. CI biochemical reaction model fit
% the flux for CITS reaction in this model was devolped in Audi & Dash lab
% previously(Xiao et al 2018)
% Model fitted to data from Bazil et al. 2010
% The units are as follows Concentration M, Flux mmol/nin, pH= 8.1, Volume =1 ml, Mass in Microg,5 parameters


%Driver
mpar = load('mpar.txt');

Data = load('Data.txt'); 
Flux_Data = Data(:,1); %Flux%
Ca_data = Data(:,2); %Ca 
Mg_data = Data(:,3); %Mg

%figure 1 
Ca=[0:0.01:35];
for i=1:1:length(Ca)
    Ca_tot = Ca(i)/10^6;
    Conc = [Ca_tot, 0];
    flux1(i,:)=Model(mpar,Conc);
end
h1 = figure(1)
set(h1,'Position',[10 10 300 300])
%subplot(1,2,1);
plot(Ca_data(1:6,:),Flux_Data(1:6,:),'k*',Ca,flux1,'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('Ca(\muM)')
ylabel('%Activity')
ylim([0 1.1])
text(2,1.1,'B','Fontsize',16)
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure 1 
% Mg=[5:0.01:1100];
% for i=1:1:length(Mg)
%     Mg_tot = Mg(i)/10^6;
%     Conc = [Mg_tot, 0];
%     flux2(i,:)=Model(mpar,Conc);
% end
% subplot(1,2,2);
% plot(log10(Mg_data(7:13,:)),Flux_Data(7:13,:),'k*',log10(Mg),flux2,'k','LineWidth',2.5,'MarkerSize',10.0);
% set(gcf,'color','w')
% set(gca,'FontSize',16)
% xlabel('Mg(\muM)')
% ylabel('%Activity')
% text(5,40,'C','Fontsize',16)
% box off
% 
