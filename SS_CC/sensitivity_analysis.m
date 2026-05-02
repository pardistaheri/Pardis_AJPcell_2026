clear all;
clear variables; clc; close all;
warning off; format short e;

%% Load required data and parameters
addpath('/Users/pardis/Library/CloudStorage/OneDrive-mcw.edu/Computational_Model_Pardis/Lung_mito_model_norm_HT_HS/Normoxia/SS_CC/Enzyme');
load Model_Params.mat;
load Vmaxs_norm.mat;
run Read_Condt.m;

% Load experimental data
load Exp_data_Lung.mat;  
data = data_OM;

%% Set up model settings
model_settings.verbose = false;
model_settings.normalize_errors = true;
model_settings.param_lower_bounds = 1e-9 * ones(1, p.NPar);
model_settings.param_upper_bounds = 1e6 * ones(1, p.NPar);

%% Get optimal parameter values
if exist('Vmaxs_norm', 'var')
    mpar_opt = Vmaxs_norm;
else
    mpar_opt = p.pest(1:p.NPar);
end

%% Calculate baseline error (E0)
E0 = Model_obj(mpar_opt, p, data);
fprintf('Baseline error (E0): %.6f\n', E0);

%% Part 1: Detailed parameter variation analysis with 5% increments
del_range = -0.5:0.01:0.5; % ±50% with 5% increments
n_points = length(del_range);
E_variation = zeros(length(mpar_opt), n_points);

fprintf('Running detailed parameter variation analysis...\n');
for z = 1:length(mpar_opt)
    for d = 1:n_points
        mpar_temp = mpar_opt;
        variation_factor = 1 + del_range(d);
        mpar_temp(z) = max(model_settings.param_lower_bounds(z), mpar_opt(z) * variation_factor);
        E_variation(z, d) = Model_obj(mpar_temp, p, data);
%         if E_variation(z, d)<E0
%             E_variation(z, d)=(E0-E_variation(z, d))+E0;
%         end
    end
    
    if model_settings.verbose && mod(z,5)==0
        fprintf('Completed %d/%d parameters\n', z, length(mpar_opt));
    end
end

%% Part 2: Calculate normalized sensitivity coefficients
del_sensitivity = 0.01;
S_Pj = zeros(length(mpar_opt), 1);

fprintf('Calculating normalized sensitivity coefficients...\n');
for z = 1:length(mpar_opt)
    % Positive perturbation
    mpar_plus = mpar_opt;
    mpar_plus(z) = mpar_opt(z) * (1 + del_sensitivity);
    E_plus = Model_obj(mpar_plus, p, data);
    
    % Negative perturbation
    mpar_minus = mpar_opt;
    mpar_minus(z) = max(model_settings.param_lower_bounds(z), mpar_opt(z) * (1 - del_sensitivity));
    E_minus = Model_obj(mpar_minus, p, data);

    % Normalized sensitivity
   
        S_Pj(z) = ((E_plus - E_minus) / (2*del_sensitivity)) / E0;
 
end

%% Plot Part 1: E/E0 vs P/P0 with grouped subplots
clear all 
close all
clc

% del_range = -0.5:0.05:1.5; % ±50% with 5% increments
load Sensitivity_Analysis_Results.mat

mpar_opt = sensitivity_results.parameters ;
E0 = sensitivity_results.E0 ;
E_variation = sensitivity_results.E_variation ;
del_range = sensitivity_results.del_range ;
S_Pj = sensitivity_results.S_Pj ;
% sensitivity_results.parameter_names = {'PDH','CITS','ACON','IDH','AKGDH','SCAS','NDK','FH','MDH','GOT','CI','CII',...
%     'CIII','CIV','CV','AK','PYRH','GLUH','DCCS','DCCM','TCC','OME','GAE','ANT','PIC','HLeak'};


figure(1);
Position1 = [.25,.25, 3, 9];
text_size1= 10;
text_size2= 14;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;

set(figure(1),'Units','inches','Position',Position1,'PaperPosition',Position1);

% Define parameter groups
group1 = 1:10;  % Fluxes 1-10 & 16 (added AK at position 16)
group1 = [group1, 16]; % Add AK to group 1
group2 = 11:15; % Fluxes 10-15
group3 = 17:26; % Fluxes 17-26

% Subplot 1: Fluxes 1-9 & 16
subplot(3,1,1);
hold on;
colors = lines(length(group1));
for i = 1:length(group1)
    idx = group1(i);
    P_P0 =  del_range;
%     E_E0 = E_variation(idx, :) / E0;
    E_E0 = E_variation(idx, :) / min(E_variation(idx, :));
    plot(P_P0, E_E0, '-', 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'MarkerSize', 4, 'DisplayName', sensitivity_results.parameter_names{idx});
end
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
xlabel('P/P_0');
ylabel('E/E_0');
title('TCA','Fontsize',text_size2);
legend('show', 'Location', 'best', 'NumColumns', 2);
ylim([1 15])
grid off;
box off;
legend box off

% Subplot 2: Fluxes 12-16
subplot(3,1,2);
hold on;
colors = lines(length(group2));
for i = 1:length(group2)
    idx = group2(i);
    P_P0 = del_range;
%         E_E0 = E_variation(idx, :) / E0;
        E_E0 = E_variation(idx, :) / min(E_variation(idx, :));
    plot(P_P0, E_E0, '-', 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'MarkerSize', 4, 'DisplayName', sensitivity_results.parameter_names{idx});
end
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
xlabel('P/P_0');
ylabel('E/E_0');
title('ETC','Fontsize',text_size2);
legend('show', 'Location', 'best');
ylim([1 15])
grid off;
box off;
legend box off

% Subplot 3: Fluxes 17-26
subplot(3,1,3);
hold on;
colors = lines(length(group3));
for i = 1:length(group3)
    idx = group3(i);
    P_P0 = del_range;
    E_E0 = E_variation(idx, :) / min(E_variation(idx, :));
    
    plot(P_P0, E_E0, '-', 'Color', colors(i,:), 'LineWidth', 2.5, ...
         'MarkerSize', 4, 'DisplayName', sensitivity_results.parameter_names{idx});
end
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on
xlabel('P/P_0');
ylabel('E/E_0');
title('Transporters','Fontsize',text_size2);
legend('show', 'Location', 'best', 'NumColumns', 2);
ylim([1 15])
grid off;
box off;
legend box off

% Plot sensitivity coefficients in the 4th subplot
figure(2);
Position2 = [.25,.25, 6, 3];
text_size1= 12;
text_size2= 18;
linewidth= 1;
markersize=2;
MarkerSizeErr= 1;

set(figure(2),'Units','inches','Position',Position2,'PaperPosition',Position2);
set(gcf,'color','w'); set(gca,'Fontsize',text_size1,'linewidth',linewidth,'FontWeight','bold'); hold on

bar(S_Pj, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'k', 'LineWidth', 1);
set(gca, 'xtick', 1:length(mpar_opt), ...
         'xticklabel', sensitivity_results.parameter_names, ...
         'XTickLabelRotation', 45,'Fontsize',10);
set(gcf,'color','w'); set(gca,'linewidth',2,'FontWeight','bold');
ylabel('Normalized Sensitivity Coefficient');
xlabel('Parameters');
title('Control Parameter Sensitivity Coefficients','Fontsize',text_size2);
ylim ([-1.6 3.5])

grid off;
box off;


%% Display summary results
fprintf('\n=== Sensitivity Analysis Results ===\n');
fprintf('Top 5 most sensitive parameters:\n');
[~, sensitivity_rank] = sort(abs(S_Pj), 'descend');
for i = 1:5
    idx = sensitivity_rank(i);
    fprintf('%d. %s (|S_Pj| = %.3f)\n', i, sensitivity_results.parameter_names{idx}, abs(S_Pj(idx)));
end

fprintf('\nAnalysis completed. Results saved to Sensitivity_Analysis_Results.mat\n');