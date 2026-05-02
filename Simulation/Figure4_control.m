%% Single ADP figures Normoxia
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare variables & read parameters 
% X: state variables;  T: time in min;  J: fuxes;  Pk: peak values of fluxes;  HPk: half the peak values  
% c stands for cell structure; v stand for vector 
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax

% load data
% data coloumns indicate: [PM GM AM Suc SucRot]
JO2_1ADP=[data.JO2_1ADP_nmol(:,1) data.JO2_1ADP_nmol(:,2) data.JO2_1ADP_nmol(:,3) data.JO2_1ADP_nmol(:,4)]; % nmol/min/mg
O2_1ADP=[data.O2_1ADP(:,1) data.O2_1ADP(:,2) data.O2_1ADP(:,3) data.O2_1ADP(:,4)]; % microM
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP1_add=[100]*1e-6; p.ADPL1=length(ADP1_add); % uM
   options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3, ...
                    'MaxOrder', 2); % Lower order for stability
                

for i=p.ISub:1:p.NSub 
X0=ICs(p); p_tem=p; 
%%% solving ODEs and calculating state variables 
T0=0; jj=1; 
cc=1;
for ii=1:1:3 % ii=1: add substrate, ii=2-6 add 4 doses of ADP
    %%% Substrate addition 
    
    if i==4
         p_tem.ini_VTmax(p.iCI)=0*p_tem.ini_VTmax(p.iCI); % inhibit CI
    end 
    
    if ii==1
        X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
        X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
        X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
        X0(p_tem.iSUCe)=cc*SUC_index(i)*10e-3; % mM
    end
    
    %%% ADP addition
    if ii==2   
        X0(p_tem.iADPe)=X0(p.iADPe)+ADP1_add(jj);
        
    end
    
    %%% Solving ODEs
    tspan=[T0:p.tstepplot:(T0+p.time1(ii,1))];
    [T,X] = ode15s(@ODEs, tspan, X0,options, p_tem);  
        

    T0=T(end,:);        X0=X(end,:);  % redefining initial values for the next time window 
    Tc(ii,i)={T};
    Xc(ii,i)={X};  
    %%% Calculating fluxes
    for zz=1:1:length(T) %length(tspan) %length(T)
        J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
    end
    Jc(ii,i)={J(1:zz,:,i)};

end % closing ii for-loop time protocol 

% storing cell variables in vectors 
Tv(:,i)=        [Tc{1,i};   Tc{2,i};]; % min
Xv(:,:,i)=1e6*  [Xc{1,i};   Xc{2,i};]; % M
Jv(:,:,i)=1e9*  [Jc{1,i};   Jc{2,i};]./1.4; % nmol/min/mg mito
end


% plots
%%% CIV fluxes 

Position1= [.5,.5, 3, 12];
text_size1= 12;
text_size2= 18;
linewidth= 1;
MarkerSizeErr= 1;
set(figure(4),'Units','inches','Position',Position1,'PaperPosition',Position1);
cl =['r','g','b','m'];
for i =p.ISub:1:p.NSub
    subplot(4,1,1)
p.sp=1; % start of plot 
OCR(:,i)=movmean(Jv(p.sp:end,p.iCIV,i),5); % smooth the flux
plot(Tv(:,i),OCR(:,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
     
ylabel('OCR (nmol/min/mg)');
xlabel('Time (min)');
title('Control','Fontsize',18,'FontWeight','bold')
legend('P+M','G+M','Suc','Suc+Rot','Fontsize',text_size1)
xlim([0 10])
ylim([0 50])
text(.5,50,'A','Fontsize',text_size2,'FontWeight','bold')
legend box off
box off
end
for i =p.ISub:1:p.NSub
    subplot(4,1,2)

plot(Tv(:,i),JO2_1ADP(:,i),cl(i),'linewidth',linewidth,'linestyle','-.','linewidth',2);
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
ylabel('OCR (nmol/min/mg)');
xlabel('Time (min)');
xlim([0 10])
ylim([0 50])
text(.5,50,'D','Fontsize',text_size2,'FontWeight','bold')
box off

end
for i =p.ISub:1:p.NSub
    subplot(4,1,3)
p.sp=1; % start of plot 

plot(Tv(:,i),Xv(:,37,i),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
     
ylabel('O_2 conc. (\muM)');
xlabel('Time (min)');

xlim([0 10])
ylim([140 220])
text(.5,220,'G','Fontsize',text_size2,'FontWeight','bold')
box off
end
for i =p.ISub:1:p.NSub
    subplot(4,1,4)

plot(Tv(:,i),O2_1ADP(:,i),cl(i),'linewidth',linewidth,'linestyle','-.','linewidth',2);
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
ylabel('O_2 conc. (\muM)');
xlabel('Time (min)');
xlim([0 10])
ylim([150 220])
text(.5,220,'J','Fontsize',text_size2,'FontWeight','bold')
box off

end