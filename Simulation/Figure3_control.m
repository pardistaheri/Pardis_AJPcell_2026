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

Position3= [.25,.25, 3, 6];
text_size1= 14;
text_size2= 18;
linewidth= 2;
markersize=2;
MarkerSizeErr= 1;
cl =['r','g','b','m'];
set(figure(3),'Units','inches','Position',Position3,'PaperPosition',Position3)

for i=p.ISub:1:p.NSub 
    subplot(2,1,1)
if i==1
plot(Tv(:,i),dPsi(:,1),cl(i),'linestyle','-.','linewidth',linewidth); hold on
elseif i==3
plot(Tv(:,i),dPsi(:,3),cl(i),'linestyle','-.','linewidth',linewidth); hold on
end

set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
     hold on
legend('P+M','Suc','Suc+Rot','Fontsize',text_size1,'FontWeight','bold')
title('Control','Fontsize',18,'FontWeight','bold')
text(1,200,'A','Fontsize',text_size2,'FontWeight','bold')
ylabel('Measured \Delta\Psi_m (mV)','Fontsize',14,'FontWeight','bold')
xlabel('Time (min)','Fontsize',16,'FontWeight','bold');
ylim([80 200])
yticks([80 120 160 200])
xticks([5 10 15 20])
legend box off
box off;
end

p.sp=1; % start of plot 
subplot(2,1,2)
plot(Tv(30:60,1),85,'r','linewidth',2); hold on
plot(Tv(30:60,3),85,'b','linewidth',2); hold on
plot(Tv(61:end,1),Xv(61:end,p.idPsi,1)-15,'r','linewidth',2); hold on
plot(Tv(61:end,3),Xv(61:end,p.idPsi,3)-20,'b','linewidth',2); hold on
xlim([1 20])
xticks([5 10 15 20])
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold');
hold on
ylabel('Model predicted \Delta\Psi_m (mV)','Fontsize',14,'FontWeight','bold')
xlabel('Time (min)','Fontsize',16,'FontWeight','bold');
text(2,200,'D','Fontsize',text_size2,'FontWeight','bold')
yticks([80 120 160 200])
xticks([5 10 15 20])
ylim([80 200])
box off;
