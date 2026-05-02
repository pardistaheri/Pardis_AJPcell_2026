%% Sequential ADP simulation Normoxia
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


% Read data
JO2_Pkd=[data.JO2_Pk(:,1) data.JO2_Pk(:,2) data.JO2_Pk(:,3) data.JO2_Pk(:,4)]; % nM
JO2_Pk_Err=[data.JO2_Pk_Err(:,1) data.JO2_Pk_Err(:,2) data.JO2_Pk_Err(:,3) data.JO2_Pk_Err(:,4)]; % nM
JO2_nmol=[data.JO2_nmol_TC(:,1) data.JO2_nmol_TC(:,2) data.JO2_nmol_TC(:,3) data.JO2_nmol_TC(:,4)]; % nM
dPsi=[data.dPsi2(:,1) data.dPsi2(:,2) data.dPsi2(:,3)]; % nM
RCI_data=[data.RCI(:,1) data.RCI(:,2) data.RCI(:,3) data.RCI(:,4)]; % nM
RCI_Err=[data.RCI_Err(:,1) data.RCI_Err(:,2) data.RCI_Err(:,3) data.RCI_Err(:,4)]; % nM

    % substrates and ADP additions 
    PYR_index = [1 0 0 0];
    GLU_index = [0 1 0 0];
    MAL_index = [1 1 0 0];
    SUC_index = [0 0 1 1];
    ADP_add = [12.5 25 50 100] * 1e-6;
    p.ADPL = length(ADP_add);
    
    options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
    
    
    % Preallocate cell arrays
    Tc = cell(p.ADPL + 3, p.NSub); % p.Es = 3
    Xc = cell(p.ADPL + 3, p.NSub);
    Jc = cell(p.ADPL + 3, p.NSub);
    Pkc = cell(p.ADPL + 3, p.NSub);
    HPkc = cell(p.ADPL + 3, p.NSub);
    RCRc = cell(p.ADPL + 3, p.NSub);
    
    for i = p.ISub:1:p.NSub 
        X0 = ICs(p); 
        p_tem = p; 
        
        T0 = 0; %start of time
        jj = 1; %ADP addition counter
        Es = 3; % extra states (state 1, state 2, and state 4) 
        cc = 1;
        
        % Initialize Pk array for this substrate
        Pk = zeros(p.ADPL + Es, 1);
        HPk = zeros(p.ADPL + Es, 1);
        RCR = zeros(p.ADPL + Es, 1);
        
        for ii = 1:(p.ADPL + Es)
            %%% Substrate addition 
            if i == 4 % Suc+Rot, inhibition of ComplexI from the start of experiment
                p_tem.ini_VTmax(p.iCI) = 0 * p_tem.ini_VTmax(p.iCI);
            end 
            
            if ii == (Es - 1)
                X0(p_tem.iMALe) = cc * MAL_index(i) * 5e-3;
                X0(p_tem.iPYRe) = cc * PYR_index(i) * 10e-3;
                X0(p_tem.iGLUe) = cc * GLU_index(i) * 10e-3;
                X0(p_tem.iSUCe) = cc * SUC_index(i) * 7e-3;
            end
            
            %%% ADP addition
            if ii >= Es && ii ~= (p.ADPL + Es)
                X0(p_tem.iADPe) = X0(p.iADPe) + ADP_add(jj);
                jj = jj + 1;
            end
            
            %%% Solving ODEs
            tspan = T0:p.tstep:(T0 + p.time(ii,1));
            [T, X] = ode15s(@ODEs, tspan, X0, options, p_tem);  
            
                
            T0 = T(end);
            X0 = X(end,:);
            
            
            Tc{ii,i} = T;
            Xc{ii,i} = X;
            
            %%% Calculating fluxes - FIXED HERE
            n_fluxes = length(Fluxes(X(1,:), p_tem)); % Get number of fluxes
            J = zeros(length(T), n_fluxes); % Preallocate J with correct dimensions
            
            for zz = 1:length(T)
                flux_result = Fluxes(X(zz,:), p_tem);
                J(zz,:) = flux_result; 
            end
            
            Jc{ii,i} = J;
            
            s = 14; % CIV flux- OCR
            st = 20; 
            
            %%% finding peaks of state 3 & RCI
            if ii < Es 
                Pk(ii,i) = mean(J(st:end-20, s));
                HPk(ii,i) = mean(J(st:end-20, s));
                RCR(ii,i) = 1;
                XPk(ii,:,i)=mean(Xc{ii,i});
                Xmin(ii,:,i)=mean(Xc{ii,i});
            else
                Pk(ii,i) = max(J(:, s));
                HPk(ii,i) = 0.5 * Pk(ii,i);
                XPk(ii,:,i)=max(Xc{ii,i});
                Xmin(ii,:,i)=min(Xc{ii,i});
                if ii > 2 % Make sure we have a valid denominator
                    RCR(ii,i) = Pk(ii,i)*1.35 / Pk(2,i);
                else
                    RCR(ii,i) = 1;
                end
            end
            
            Pkc{ii,i} = Pk(ii,i);
            HPkc{ii,i} = HPk(ii,i);
            RCRc{ii,i} = RCR(ii,i);
        end
        
        % storing cell variables in vectors 
Tv(:,i)=        [Tc{1,i};   Tc{2,i};    Tc{3,i};    Tc{4,i};    Tc{5,i};    Tc{6,i};]; % min
Xv(:,:,i)=1e0*  [Xc{1,i};   Xc{2,i};    Xc{3,i};    Xc{4,i};    Xc{5,i};    Xc{6,i};]; % M
Jv(:,:,i)=1e9*  [Jc{1,i};   Jc{2,i};    Jc{3,i};    Jc{4,i};    Jc{5,i};    Jc{6,i};]./1.4; % nmol/min/mg mito
JPkv(:,:,i)=1e9*[Pkc{1,i}; Pkc{2,i};  Pkc{3,i};  Pkc{4,i};  Pkc{5,i};  Pkc{6,i};]./1.4; % nmol/min/mg mito
RCRv(:,:,i)=[RCRc{1,i}; RCRc{2,i}; RCRc{3,i}; RCRc{4,i}; RCRc{5,i}; RCRc{6,i};]; %unitless
JO(:,i)=Jv(:,25,i);
    end

%%%%%%%%%%%Figure 4
Position4= 3*[.25,.25, 5.5, 5];
text_size1= 10;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];

set(figure(4),'Units','inches','Position',Position4,'PaperPosition',Position4)
titles = ["PDH","CITS","ACON","IDH","AKGDH","SCAS","NDK","FH",...
    "MDH","GOT","CI","CII","CIII","CIV","CV","AK","PYRH","GLUH","DCCS",...
    "DCCM","TCC","OME","GAE","ANT","PIC","HLEAK"];
p.sp=1;
for i =p.ISub:1:p.NSub
   for j= 1:1:26
     subplot(7,4,j)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,j,i),cl(i),'linewidth',linewidth)
%      xlabel('Time (sec)');   

title(titles(j));
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     hold on
     xlim([3 20])
     ylabel('J(nmol/min/mg)','Fontsize',10,'FontWeight','bold')
     xlabel('Time (min)','Fontsize',10,'FontWeight','bold');
% ax=gca; ax.XLim=[0 inf];
box off;  
   end
end

%%%%%%%%%%% Figure 5
set(figure(5),'Units','inches','Position',Position4,'PaperPosition',Position4);


titles = ["ADPe","ATPe","Pie","PYRe","MALe","CITe",'AKGe',"SUCe",...
    "GLUe","ASPe","He","ADPm","ATPm","AMPm","GDPm","GTPm","NADm","NADHm","UQm",...
    "UQH2m","CytCox","CytCred","Pim","GLUm","ASPm","PYRm",...
    "OXAm","CITm","ICITm","AKGm","SCOAm","SUCm","FUMm",...
    "MALm","COAm","ACOAm","O2m","dPsim","Hm"];


for i =p.ISub:1:p.NSub
   for j= 1:1:6
      subplot(6,6,j)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,j,i),cl(i),'linewidth',linewidth); hold on
    title(titles(j)); 
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     hold on
     ylabel('Conc.(M)','Fontsize',10,'FontWeight','bold')
     xlabel('Time (min)','Fontsize',10,'FontWeight','bold')
     xlim([3 20])
     box off; 
    for k= 8:1:37
        subplot(6,6,k-1)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,k,i),cl(i),'linewidth',linewidth); hold on
     title(titles(k));
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     hold on
     ylabel('Conc.(M)','Fontsize',10,'FontWeight','bold')
     xlabel('Time (min)','Fontsize',10,'FontWeight','bold')
     xlim([2 20])
% ax=gca; ax.XLim=[0 inf];
box off; 
   
end
   end
end

%%%%%%%%%%% Figure 7
Position7= [.25,.25, 3, 15];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
set(figure(7),'Units','inches','Position',Position7,'PaperPosition',Position7);
 
for i =3:1:p.NSub
    subplot(5,1,1)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,32,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('SUC_m(\muM)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     title('Control','Fontsize',18,'FontWeight','bold')
     legend('Suc','Suc+Rot','Fontsize',text_size1,'FontWeight','bold')
     text(4,30000,'A','Fontsize',text_size2,'FontWeight','bold')
     ylim([10000 30000]) 
     xlim([3 20])
     legend box off
     box off
     subplot(5,1,2)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,33,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('MAL_m(\muM)','Fontsize',12,'FontWeight','bold'); 
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold');
     text(4,400,'D','Fontsize',text_size2,'FontWeight','bold')
     ylim([100 400])
     xlim([3 20])
     box off;
     
     subplot(5,1,3)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,27,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('OXA_m(\muM)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     text(4,0.001,'G','Fontsize',text_size2,'FontWeight','bold')
     ylim([0 0.001])
     xlim([3 20])
     box off;
     
     subplot(5,1,4)
     plot(Tv(p.sp:end,i),Xv(p.sp:end,23,i).*10^6,cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('Pi_m(\muM)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     text(4,1.5*10^4,'J','Fontsize',text_size2,'FontWeight','bold')
     ylim([1 1.5]*10^4)
     xlim([3 20])
     box off;
     
     subplot(5,1,5)
     plot(Tv(30:60,3),85,'b','linewidth',2); hold on
     plot(Tv(30:60,4),85,'m','linewidth',2); hold on
     
     plot(Tv(61:end,3),Xv(61:end,p.idPsi,3)-20,'b','linewidth',2); hold on
     plot(Tv(61:end,4),Xv(61:end,p.idPsi,4)-20,'m','linewidth',2); hold on

%      plot(Tv(p.sp:end,i),Xv(p.sp:end,38,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('\Delta\Psi_m(mV)','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')  
     text(4,200,'M','Fontsize',text_size2,'FontWeight','bold')
     ylim([100 200])
     xlim([3 20])
     box off;
end
 
%%%%%%%%%%% Figure 8
Position8= [.25,.25, 3, 15];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
set(figure(8),'Units','inches','Position',Position8,'PaperPosition',Position8);
 
for i =3:1:p.NSub
    subplot(5,1,1)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,19,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_D_C_C_S','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     title('Control','Fontsize',18,'FontWeight','bold')
     legend('Suc','Suc+Rot','Fontsize',text_size1,'FontWeight','bold')
     text(4,60,'A','Fontsize',text_size2,'FontWeight','bold')
     ylim([-20 60])
     xlim([3 20]) 
     legend box off
     box off
     subplot(5,1,2)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,20,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_D_C_C_M','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     text(4,0,'D','Fontsize',text_size2,'FontWeight','bold')
     ylim([-50 0])
     xlim([3 20])
     box off;
     subplot(5,1,3)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,12,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_S_D_H','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     text(4,50,'G','Fontsize',text_size2,'FontWeight','bold')
     ylim([0 50])
     xlim([3 20])
     box off;
     subplot(5,1,4)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,25,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_P_I_C','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')
     text(4,70,'J','Fontsize',text_size2,'FontWeight','bold')
     ylim([-10 70])
     xlim([3 20])
     box off;
     subplot(5,1,5)
     plot(Tv(p.sp:end,i),Jv(p.sp:end,14,i),cl(i),'linewidth',2); hold on
     set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth);
     ylabel('J_O_2','Fontsize',12,'FontWeight','bold');   
     xlabel('Time(min)','Fontsize',12,'FontWeight','bold')  
     text(4,50,'M','Fontsize',text_size2,'FontWeight','bold')
     ylim([0 50])
     xlim([3 20])
     box off;
end
