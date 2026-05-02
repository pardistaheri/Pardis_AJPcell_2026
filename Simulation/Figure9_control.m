%% clear workspace 
clear all
close all
clc
warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read parameters
run Read_Condt
load Vmaxs_norm.mat
p.pest(1:p.NPar)=Vmaxs_norm; p.pest; % estimated Vmax
% experimental condition
p.time=[2, 1, 1, 5,1]';
p.NSub=2; % number of substrates  and last substrate
% substrates and ADP additions 
PYR_index=[1 0 0 0];
GLU_index=[0 1 0 0];
MAL_index=[1 1 0 0];
SUC_index=[0 0 1 1];
ADP_add=[0 100]*1e-6; p.ADPL=length(ADP_add); % uM    
options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
CIF =[1:-0.01:0]; % CI factor
counter1=0; p.Es=3; counter4=20;
counter5=1; % to store T, X and J with different lengths
T0ii=zeros(p.ADPL+p.Es,5);       X0ii=zeros(p.ADPL+p.Es,5);

% Pre-allocate arrays for storing results
JPk_O2 = zeros(length(CIF), p.NSub);
JPk_ANT = zeros(length(CIF), p.NSub);
JPk_CII = zeros(length(CIF), p.NSub);
JPk_CIII = zeros(length(CIF), p.NSub);
JPk_CIV = zeros(length(CIF), p.NSub);
JPk_CV = zeros(length(CIF), p.NSub);

for j=1:length(CIF)
    CIF(j);
    for i=p.ISub:1:p.NSub 
        X0=ICs(p); p_tem=p; 
        p_tem.ini_VTmax(p.iCI)=CIF(j)*p_tem.ini_VTmax(p.iCI);
        
        %%% solving ODEs and calculating state variables 
        T0=0; jj=1; p.Es=3; % extra states 
        cc=1;
        
        for ii=1:1:p.ADPL+p.Es 
            %%% Substrate addition       
%             if i == 4 % Suc+Rot, inhibition of ComplexI from the start of experiment
%                 % Only apply Rotenone inhibition, don't set CI to 0
%                 % p_tem.ini_VTmax(p.iCI) = 0 * p_tem.ini_VTmax(p.iCI);
%                 % COMMENT OUT OR REMOVE the line above to allow full CI range
%             end
            
            if ii==p.Es-1
                X0(p_tem.iMALe)=cc*MAL_index(i)*5e-3; % mM
                X0(p_tem.iPYRe)=cc*PYR_index(i)*10e-3; % mM
                X0(p_tem.iGLUe)=cc*GLU_index(i)*10e-3; % mM
                X0(p_tem.iSUCe)=cc*SUC_index(i)*7e-3; % mM
            end
            
            %%% ADP addition
            if ii>=p.Es && ii~=p.ADPL+p.Es   
                X0(p_tem.iADPe)=X0(p.iADPe)+ADP_add(jj);
                jj=jj+1;
            end
            
            %%% Solving ODEs
            tspan=[T0:p.tstep:(T0+p.time(ii,1))];
            [T,X] = ode15s(@ODEs, tspan, X0, options, p_tem);      
            T0=T(end,:);       X0=X(end,:);  
            TL(ii,i,j)=length(T);
            Tc(ii,i,j)={T};      Xc(ii,i,j)={X}; 
            
            %%% Calculating fluxes
            for zz=1:1:length(T)
                J(zz,:,i)=Fluxes(X(zz,:),p_tem);        
            end
            Jc(ii,i,j)={J(1:zz,:,i)};

            %%% finding peaks of state 3 for ODEs
            s=25; % CIV flux- OCR
            st=20; 
            
            %%% finding peaks and length of state 3
            if ii<p.Es 
                Pk(ii,i,j)=mean(J(st:end-20,s,i));  
                HPk(ii,i,j)=mean(J(st:end-20,s,i));    
                XPk(ii,i,j)=mean(X(:,i));
                RCR(ii,i,j) = 1;
            elseif ii>=p.Es || ii~=p.ADPL+p.Es
                Pk(ii,i,j)=max(J(1:end,s,i)); 
                HPk(ii,i,j)=.5*max(J(1:end,s,i));
                XPk(ii,i,j)=max(X(:,i));
                RCR(ii,i,j)=Pk(ii,i,j)/Pk(2,i,j);
            end
   
            if ii==4
                Xmean_NADH(j,i)=mean(X(10:20,18));
                Xmean_UQH2(j,i)=mean(X(10:20,20));     
                Xmean_Cytr(j,i)=mean(X(10:20,22));       
                Xmin_dpsi(j,i)= min(X(10:end,38));
                Xmin_ATP(j,i)= min(X(10:end,13));
                Xmax_ATP(j,i)= max(X(:,13));
                JPk_ANT(j,i)= 1e9*max(J(:,25))./1.4;  % ANT activity
                JPk_CII(j,i)= 1e9*max(J(:,12))./1.4;  % ANT activity
                JPk_CIII(j,i)= 1e9*max(J(:,13))./1.4;  % ANT activity
                JPk_CIV(j,i)= 1e9*max(J(:,14))./1.4;  % ANT activity
                JPk_CV(j,i)= 1e9*max(J(:,15))./1.4;  % ANT activity
                JPk_O2(j,i)= 1e9*max(J(:,14))./1.4;   % JO2
            end

            JPkc(ii,i)={Pk(ii,i)};
            HPkc(ii,i)={HPk(ii,i)};
            RCIc(ii,i)={RCR(ii,i)};
        end % for ii

        % storing cell variables in vectors 
        Tv(:,i,j)= [Tc{2,i,j}; Tc{3,i,j}; Tc{4,i,j}; Tc{5,i,j};];
        Xv(:,:,i)=1e0* [Xc{2,i,j}; Xc{3,i,j}; Xc{4,i,j}; Xc{5,i,j};]; % M
        Jv(:,:,i)=1e9* [Jc{2,i,j}; Jc{3,i,j}; Jc{4,i,j}; Jc{5,i,j};]./1.4;
    end % i substrate for-loop
end % j for leak factor

Position2= [.25,.25, 3, 15];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
CI_activity =CIF*100;
set(figure(14),'Units','inches','Position',Position2,'PaperPosition',Position2)

for i=p.ISub:1:p.NSub 
subplot(5,1,1)
plot(CI_activity,(Xmean_NADH(:,i)./3e-03),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');   
ylabel('NADH ratio','Fontsize',12,'FontWeight','bold');
xlabel('CI Activity (%)','Fontsize',12,'FontWeight','bold');
title('Control','Fontsize',18,'FontWeight','bold')
legend('P+M','G+M','Suc','Suc+Rot','Fontsize',text_size1,'FontWeight','bold')
text(95,1,'A','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
ylim([0 1]);
legend box off
box off;
end 

for i=p.ISub:1:p.NSub 
subplot(5,1,2)
plot(CI_activity,(Xmean_UQH2(:,i)./1.5e-03),cl(i),'linewidth',2); hold on  
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');  
ylabel('UQH_2 ratio','Fontsize',12,'FontWeight','bold');
xlabel('CI Activity','Fontsize',12,'FontWeight','bold');
ylim([0 0.6])
text(95,0.6,'D','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
box off;
end 

for i=p.ISub:1:p.NSub
subplot(5,1,3)
plot(CI_activity,(Xmean_Cytr(:,i)./3e-03),cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
set(gca, 'XDir', 'reverse');  
ylabel('CytC_r_e_d ratio','Fontsize',12,'FontWeight','bold');
xlabel('CI Activity','Fontsize',12,'FontWeight','bold');
ylim([0 0.6])
text(5,0.6,'G','Fontsize',text_size2,'FontWeight','bold')
xtickformat('percentage'); 
box off;
end 

for i=p.ISub:1:p.NSub
subplot(5,1,4)
if i==1||i==2
plot(CI_activity,Xmin_dpsi(:,i)-15,cl(i),'linewidth',2); hold on
else 
plot(CI_activity,Xmin_dpsi(:,i)-20,cl(i),'linewidth',2); hold on
end
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on 
set(gca, 'XDir', 'reverse');  
ylabel('\Delta\Psi_m (mV)','Fontsize',12,'FontWeight','bold');
xlabel('CI Activity','Fontsize',12,'FontWeight','bold');
xtickformat('percentage'); 
ylim([120 160])
text(95,160,'J','Fontsize',text_size2,'FontWeight','bold')
box off;
end

for i=p.ISub:1:p.NSub
subplot(5,1,5)
plot(CI_activity,(Xmax_ATP(:,i)-Xmin_ATP(:,i))*10^3,cl(i),'linewidth',2); hold on
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on 
set(gca, 'XDir', 'reverse');  
ylabel('ATP_m (mM)','Fontsize',12,'FontWeight','bold');
xlabel('CI Activity','Fontsize',12,'FontWeight','bold');
xtickformat('percentage'); 
ylim([0 3.5])
text(95,3,'M','Fontsize',text_size2,'FontWeight','bold')
box off;
end

