% Load experimental data
Data = load('Data.txt');
Flux_Data = Data(:,1); %Flux(mM/min)
A = Data(:,2); %C_Sucm (uM)
% B = Data(:,3)/10^3; %C_FADmm (uM)
% C = Data(:,4)/10^3; %C_Fumm (uM)
% D = Data(:,5)/10^3; %C_FADH2mm (uM)

mpar = load('mpar.txt');
options = [];       % ODE options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 1
SUC_tot=[0:0.01:2.1]; %
FAD_tot=[0.5 0.67 1 2 5]*10^-3; %
FUM_tot=[0 0 0 0 0]; %
FADH2_tot=[0 0 0 0 0]; %
for i=1:1:length(FAD_tot)
    for j=1:1:length(SUC_tot)
        A_tot = SUC_tot(j)*10^-3;
        B_tot = FAD_tot(i);
        C_tot = FUM_tot(i);
        D_tot = FADH2_tot(i);
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux2(i,j)=flux(mpar,Conc)*10^3;
    end
end
h1 = figure(1)
set(h1,'Position',[10 10 600 300])
subplot(1,2,1);
plot( A(20:23,:),Flux_Data(20:23,:),'k*',A(24:27,:),Flux_Data(24:27,:),'ko',A(28:31,:),Flux_Data(28:31,:),'kd',A(32:35,:),Flux_Data(32:35,:),'k<', A(36:39,:),Flux_Data(36:39,:),'k+',...
    SUC_tot,flux2(1,:),'k',SUC_tot,flux2(2,:),'k',SUC_tot,flux2(3,:),'k',SUC_tot,flux2(4,:),'k', SUC_tot,flux2(5,:),'k','LineWidth',1.5,'MarkerSize',8.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('SUC (mM)')
ylabel('CII/SDH Forward Flux (\mumol/min)')
legend('0.5 mM UQ', '0.67 mM UQ','1 mM UQ', '2 mM UQ','5 mM UQ','Location','northwest')
text(0.1,2.5,'A','Fontsize',18)
legend boxoff
box off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = load('Data_R.txt');
Flux_Data = Data(:,1); %Flux (nmol/min)
D = Data(:,3); %C_Sucm (uM)

mpar = load('mpar_R.txt');
options = [];       % ODE options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 2
FADH2_tot=[0:0.01:20]; %
FUM_tot=10/10^3; %
for i=1:1:length(FADH2_tot)
    A_tot =0;
    B_tot = 0;
    C_tot = FUM_tot;
    D_tot = FADH2_tot(i)/10^6;
    Conc=[A_tot, B_tot, C_tot, D_tot];
    flux2(i,:)=flux_R(mpar,Conc)*10^9;
end
subplot(1,2,2);
plot(FADH2_tot,flux2,'k', D(1:4,:),Flux_Data(1:4,:),'k*', 'LineWidth',2.0,'MarkerSize',10.0);
set(gcf,'color','w')
set(gca,'FontSize',14,'LineWidth',2.0)
xlabel('UQH_2 (\muM)')
ylabel('CII/SDH Reverse Flux (nmol/min)')
text(1,1,'B','Fontsize',18)
box off
