% Load experimental data
Data = load('Data_R.txt');
Flux_Data = Data(:,1); %Flux(mM/min)
D = Data(:,3); %C_Sucm (uM)
% B = Data(:,3)/10^3; %C_FADmm (uM)
% C = Data(:,4)/10^3; %C_Fumm (uM)
% D = Data(:,5)/10^3; %C_FADH2mm (uM)

mpar = load('mpar_R.txt');
options = [];       % ODE options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
FADH2_tot=[0:0.01:20]; %
FUM_tot=10/10^3; %
for i=1:1:length(FADH2_tot)
    A_tot =0;
    B_tot = 0;
    C_tot = FUM_tot;
    D_tot = FADH2_tot(i)/10^6;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux1(i,:)=flux_R(mpar,Conc)*10^9;
end
plot(FADH2_tot,flux1,'k', D(1:4,:),Flux_Data(1:4,:),'k*', 'LineWidth',2.0,'MarkerSize',10.0);set(gcf,'color','w');
set(gca,'FontSize',12)
xlabel('FADH2 (\muM)')
ylabel('Flux (nmol/min)')
box off
