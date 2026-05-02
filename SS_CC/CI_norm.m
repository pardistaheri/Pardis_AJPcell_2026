%% Main script for parameter correlation analysis
clear all; close all; clc;

% Load necessary data and parameters
addpath('/Users/pardis/Library/CloudStorage/OneDrive-mcw.edu/Computational_Model_Pardis/Lung_mito_model_norm_HT_HS/Normoxia/SS_CC/Enzyme');
load Model_Params;
load Exp_data_Lung;  
data = data_OM;
JO2_data=data.CI; % nM
% Load best-fit parameters (replace with your actual best-fit parameters)
load Vmaxs_norm; % This should contain the best-fit parameter values
p.pest = Vmaxs_norm; % Adjust based on your actual variable name

% Initialize parameters
p.org = 1; % organ 
p.Algt = 1; % algorithm 
p.ISub = 1; % initial substrate 
p.NSub = 3; % number of substrates and last substrate
p.NPar = 26; % number of parameters
p.NOde = 39; % number of ODEs
p.close_system = 1; % 1 simulates a closed system like Oroboros, and 0 simulates an open system like PTI
p.Q10_corr = p.Q10_con.^((p.Tem-p.Tem_Stnd)/10);
p.tstep = 1/30;   % time step min
p.tstepplot = 1/30;   % time step min
p.beta = 0.35;
p.Ve = 2e-3; % buffer_Vol/mito_Vol experimet:.2mg/mL
p.time = 1*[2, 3, 2, 3, 4, 4, 2]'; % Seq ADP time

% parameter values
p0 = p.pest(1:p.NPar); % Best-fit parameters

% Compute JO2 flux
[JO2_nominal, ~] = compute_JO2(p0, p);
Nobs = numel(JO2_nominal); % Total number of observations
NPar = length(p0); % Number of parameters

% Initialize Jacobian matrix
JM = zeros(Nobs, NPar);

% Compute Jacobian using central difference
fprintf('Computing Jacobian matrix...\n');
for j = 1:NPar
    fprintf('Processing parameter %d of %d\n', j, NPar);
    
    % 1% of parameter value
    dp = 0.01 * p0(j);
    
    % Positive perturbation
    p_plus = p0;
    p_plus(j) = p0(j) + dp;
    [JO2_plus, ~] = compute_JO2(p_plus, p);
    
    % Negative perturbation
    p_minus = p0;
    p_minus(j) = p0(j) - dp;
    [JO2_minus, ~] = compute_JO2(p_minus, p);
    
    % Central difference approximation
    JM(:, j) = (JO2_plus - JO2_minus) / (2 * dp); 
end


%% Calculate confidence intervals using nlparci
fprintf('Computing parameter confidence intervals...\n');

    residuals = (JO2_nominal - JO2_data)/100; save rsiduals;

% Calculate confidence intervals (95% by default)
alpha = 0.05; % 95% confidence level

% Use regularized Jacobian
CI = nlparci(p0, residuals, 'jacobian', JM, 'alpha', alpha); save CI;

% CI = nlparci(p0, residuals, 'jacobian', JM, 'alpha', alpha); save CI;

% Display confidence interval results
fprintf('\nParameter Confidence Intervals (%.1f%%):\n', (1-alpha)*100);
fprintf('Parameter #\tLower Bound\tEstimate\tUpper Bound\tRelative Uncertainty (%%)\n');
fprintf('--------------------------------------------------------------------------\n');
% 
for i = 1:NPar
    lower_bound = CI(i, 1);
    upper_bound = CI(i, 2);
    estimate = p0(i);
    rel_uncertainty = 100 * (upper_bound - lower_bound) / (2 * abs(estimate));
    
    fprintf('%d\t\t%.4e\t%.4e\t%.4e\t%.2f\n', ...
            i, lower_bound, estimate, upper_bound, rel_uncertainty);
end

% Identify parameters with high uncertainty (e.g., > 50%)
high_uncertainty_threshold = 50;
high_uncertainty_params = [];
fprintf('\nParameters with high uncertainty (>%.0f%%):\n', high_uncertainty_threshold);
for i = 1:NPar
    rel_uncertainty = 100 * (CI(i, 2) - CI(i, 1)) / (2 * abs(p0(i)));
    if rel_uncertainty > high_uncertainty_threshold
        high_uncertainty_params = [high_uncertainty_params; i rel_uncertainty];
        fprintf('Parameter %d: %.2f%% uncertainty\n', i, rel_uncertainty);
    end
end

% Plot confidence intervals
xvalues = {'\bfPDH','\bfCITS','\bfACON','\bfIDH','\bfAKGDH','\bfSCAS','\bfNDK','\bfFH','\bfMDH',...
    '\bfGOT','\bfCI','\bfCII','\bfCIII','\bfCIV','\bfCV','\bfAK','\bfPYRH','\bfGLUH','\bfDCCS',...
    '\bfDCCM','\bfTCC','\bfOME','\bfGAE','\bfANT','\bfPIC','\bfHLeak'};

figure(31);
set(figure(31), 'Units', 'inches', 'Position', [1, 1, 12, 5]);
errorbar(1:NPar, p0, p0 - CI(:,1), CI(:,2) - p0, 'o', 'LineWidth', 1.5, 'MarkerSize', 8);
% xlabel(xvalues);
set(gcf,'color','w'); set(gca,'Fontsize',16,'linewidth',2,'FontWeight','bold'); hold on
ylabel('Parameter Value');
title('Parameter Estimates with 95% Confidence Intervals','FontSize', 20);
grid off;
box off;

% Add parameter names if available
if exist('xvalues', 'var')
    set(gca, 'XTick', 1:NPar, 'XTickLabel', strrep(xvalues, '\bf',''));
    xtickangle(45);
end



%% Helper function to compute JO2 flux 
function [JO2_vector, times] = compute_JO2(pest, p)
    % Assign parameters
    p.pest(1:p.NPar) = pest;
    
    % Initialize outputs
    JO2_vector = [];
    times = [];
    OCR = [];
    
    % Substrate indices
    PYR_index = [1 0 0 0];
    GLU_index = [0 1 0 0];
    MAL_index = [1 1 0 0];
    SUC_index = [0 0 1 1];
    ADP_add = [12.5 25 50 100] * 1e-6;
    p.ADPL = length(ADP_add);
    
    % ODE solver options
    options = odeset('NonNegative', 1:p.NOde, ...
                    'RelTol', 1e-6, ...
                    'AbsTol', 1e-8, ...
                    'MaxStep', 0.1, ...
                    'InitialStep', 1e-3);
    
    % Loop over substrates
    for i = p.ISub:1:p.NSub
%         

        X0 = ICs(p); 
        p_tem = p; 
        
        T0 = 0; 
        jj = 1; 
        Es = 3; % extra states 
        cc = 1;
        
        % Initialize for this substrate
        JO2_substrate = [];
        times_substrate = [];
        
        for ii = 1:(p.ADPL + Es)
            % Substrate addition 
            if i == 4
                p_tem.ini_VTmax(p.iCI) = 0 * p_tem.ini_VTmax(p.iCI);
            end 
            
            if ii == (Es - 1)
                X0(p_tem.iMALe) = cc * MAL_index(i) * 5e-3;
                X0(p_tem.iPYRe) = cc * PYR_index(i) * 10e-3;
                X0(p_tem.iGLUe) = cc * GLU_index(i) * 10e-3;
                X0(p_tem.iSUCe) = cc * SUC_index(i) * 7e-3;
            end
            
            % ADP addition
            if ii >= Es && ii ~= (p.ADPL + Es)
                X0(p_tem.iADPe) = X0(p.iADPe) + ADP_add(jj);
                jj = jj + 1;
            end
            
            % Solve ODEs
            tspan = T0:p.tstep:(T0 + p.time(ii,1));
            [T, X] = ode15s(@ODEs, tspan, X0, options, p_tem);  
            
            T0 = T(end);        
            X0 = X(end,:);
            
            % Calculate fluxes
            n_fluxes = 26; % Number of fluxes from Fluxes.m
            J = zeros(length(T), n_fluxes);
            
            for zz = 1:length(T)
                flux_result = Fluxes(X(zz,:), p_tem);
                J(zz,:) = flux_result;
            end
            
            % Extract JO2 flux (index 14)
            JO2_current = 1e9*movmean(J(:, 14),5);
            JO2_substrate = [JO2_substrate; JO2_current]; 
            times_substrate = [times_substrate; T];
        end
        
        % Store for this substrate
        JO2_vector = [JO2_vector; JO2_substrate]; save JO2_vector;
        times = [times; times_substrate];
        
    end
end