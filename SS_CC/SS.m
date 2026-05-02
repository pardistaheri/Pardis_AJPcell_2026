clear variables; clc; close all;
warning off; format short e;

%% Load required data and parameters
addpath('/Users/pardis/Desktop/Pardis_Lung_mito_model_norm_HT_HS/Normoxia/SS_CC/Enzyme');
load Model_Params.mat;
load Vmaxs_norm.mat;
run Read_Condt.m;

% Load experimental data
load Exp_data_Lung.mat;  
data = data_OM;

%% Set up model settings
model_settings.verbose = false;
model_settings.normalize_errors = true;
model_settings.param_lower_bounds = 1e-9 * ones(1, p.NPar); % Conservative lower bounds
model_settings.param_upper_bounds = 1e6 * ones(1, p.NPar);  % Conservative upper bounds
model_settings.error_weights = struct('OCR', 0.4, 'MMP', 0.3, 'ATP', 0.2, 'redox', 0.1);

%% Get optimal parameter values
% Assuming optimal parameters are stored in Vmaxs_norm or need to be extracted
if exist('Vmaxs_norm', 'var')
    mpar_opt = Vmaxs_norm; % Use stored optimal parameters
else
    % If parameters need to be extracted from p structure
    mpar_opt = p.pest(1:p.NPar);
end

%% Calculate baseline error (E0)
E0 = Model_obj(mpar_opt, p, data);
fprintf('Baseline error (E0): %.6f\n', E0);

%% Part 1: Impact of ±50% parameter variation on objective function
del_Pfact = 0.5; % 50% variation
E_plus50 = zeros(length(mpar_opt), 1);
E_minus50 = zeros(length(mpar_opt), 1);

fprintf('Running ±50%% parameter variation analysis...\n');
for z = 1:length(mpar_opt)
    % +50% variation
    mpar_temp = mpar_opt;
    mpar_temp(z) = mpar_opt(z) * (1 + del_Pfact);
    E_plus50(z) = Model_obj(mpar_temp, p, data);
    
    % -50% variation
    mpar_temp = mpar_opt;
    mpar_temp(z) = max(model_settings.param_lower_bounds(z), mpar_opt(z) * (1 - del_Pfact));
    E_minus50(z) = Model_obj(mpar_temp, p, data);
    
    if model_settings.verbose
        fprintf('Parameter %d: E(+50%%) = %.6f, E(-50%%) = %.6f\n', z, E_plus50(z), E_minus50(z));
    end
end

%% Part 2: Calculate normalized sensitivity coefficients using central difference
del_sensitivity = 0.01; % 1% change for central difference
S_Pj = zeros(length(mpar_opt), 1);

fprintf('Calculating normalized sensitivity coefficients...\n');
for z = 1:length(mpar_opt)
    % Positive  (+1%)
    mpar_plus = mpar_opt;
    mpar_plus(z) = mpar_opt(z) * (1 + del_sensitivity);
    E_plus = Model_obj(mpar_plus, p, data);
    
    % Negative (-1%)
    mpar_minus = mpar_opt;
    mpar_minus(z) = max(model_settings.param_lower_bounds(z), mpar_opt(z) * (1 - del_sensitivity));
    E_minus = Model_obj(mpar_minus, p, data);


% Should be:
if E0 > 0
    S_Pj(z) = ((E_plus - E_minus) / (2*del_sensitivity)) / E0;
else
    S_Pj(z) = 0;
end
       
    if model_settings.verbose
        fprintf('Parameter %d: S_Pj = %.6f\n', z, S_Pj(z));
    end
end

%% Save sensitivity analysis results
sensitivity_results.parameters = mpar_opt;
sensitivity_results.E0 = E0;
sensitivity_results.E_plus50 = E_plus50;
sensitivity_results.E_minus50 = E_minus50;
sensitivity_results.S_Pj = S_Pj;
sensitivity_results.parameter_names = {'PDH','CITS','ACON','IDH','AKGDH','SCAS','NDK','FH','MDH','GOT','CI','CII',...
    'CIII','CIV','CV','AK','PYRH','GLUH','DCCS','DCCM','TCC','OME','GAE','ANT','PIC','HLeak'};

save('Sensitivity_Analysis_Results.mat', 'sensitivity_results');

%% Plot results - Normalized objective function vs parameter variation

figure(1);
Position1= [.25,.25, 21, 12];
text_size1= 10;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;

set(figure(1),'Units','inches','Position',Position1,'PaperPosition',Position1)

for i = 1:length(mpar_opt)
    subplot(7, 4, i); % Adjusted for 26 parameters 
    
    % Plot ±50% variation points
    plot(0.5, E_minus50(i)/E0, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    plot(1.0, 1.0, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5); % Optimal point
    plot(1.5, E_plus50(i)/E0, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    
    % Add connecting lines
    plot([0.5, 1.0, 1.5], [E_minus50(i)/E0, 1.0, E_plus50(i)/E0], 'k-', 'LineWidth', 1);
    set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
    title(sensitivity_results.parameter_names{i}, 'FontSize', 15);
    xlabel('P/P_0');
    ylabel('E/E_0');
%     set(gca, 'FontSize', 7, 'LineWidth', 1);
    box off;
    grid off;
end

%% Plot normalized sensitivity coefficients
figure(2);
Position2= [.25,.25, 12, 5];
text_size1= 12;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;

set(figure(2),'Units','inches','Position',Position2,'PaperPosition',Position2)

bar(S_Pj, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'k', 'LineWidth', 1);
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
ylabel('Normalized Sensitivity Coefficient (S_{P_j})');
xlabel('Parameters');
title('Normalized Parameter Sensitivity Coefficients for Normoxia','FontSize', 18);
% ylim([-0.2 9])
set(gca, 'xtick', 1:length(mpar_opt), ...
         'xticklabel', sensitivity_results.parameter_names, ...
         'LineWidth', 1, ...
         'XTickLabelRotation', 45);
box off;
grid off;

%% Display results in table format
fprintf('\n=== Sensitivity Analysis Results ===\n');
fprintf('Parameter\t\tS_Pj\t\tE(+50%%)/E0\tE(-50%%)/E0\n');
fprintf('---------\t\t----\t\t----------\t----------\n');

for i = 1:length(mpar_opt)
    fprintf('%-10s\t\t%6.3f\t\t%8.3f\t%8.3f\n', ...
            sensitivity_results.parameter_names{i}, ...
            S_Pj(i), ...
            E_plus50(i)/E0, ...
            E_minus50(i)/E0);
end

%% Additional analysis: Identify most sensitive parameters
[~, sensitivity_rank] = sort(abs(S_Pj), 'descend');
fprintf('\nTop 5 most sensitive parameters:\n');
for i = 1:5
    idx = sensitivity_rank(i);
    fprintf('%d. %s (|S_Pj| = %.3f)\n', i, sensitivity_results.parameter_names{idx}, abs(S_Pj(idx)));
end

fprintf('\nSensitivity analysis completed. Results saved to Sensitivity_Analysis_Results.mat\n');