%% Clear worksapce
close all; 
clear all; 
clc;
warning off
format short e 
tic % measure elapsed time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters, data and variables 
run Read_Condt   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting initial conditions and boundaries 
load Params_int.mat
x0=Vmaxs_norm;  % initial Vmax

%%% Set the fmincon optimizer settings for parameter estimation
options = optimoptions('fmincon');
options = optimoptions(options,...
    'Algorithm', 'interior-point', ... % Explicitly choose a robust algorithm
    'SpecifyObjectiveGradient', false, ... % Set to true if you provide a gradient
    'Display', 'iter',...
    'MaxIterations', 1000, ...          % Drastically increase iterations
    'MaxFunctionEvaluations', 3000, ... % Increase function evaluations
    'StepTolerance', 1e-10, ...         % Replaces TolX
    'OptimalityTolerance', 1e-10, ...   % Replaces TolFun for first-order optimality
    'FunctionTolerance', 1e-10);        % Replaces TolFun for function value change

% Define Bounds
LB = x0 / 100;   % Lower bound (e.g., 1/10th of initial guess)
UB = x0 * 2;   % Upper bound (e.g., 10 times initial guess)

% multi-start optimization:
numStarts = 1;
bestParams = x0;
bestFval = inf;

for i = 1:numStarts
    % Randomize the starting point for each run to avoid local minima
    randPerturbation = 0.9 + (1.1 - 0.9) * rand(size(x0)); % Perturb between 0.9 and 1.1
    x0_perturbed = x0 .* randPerturbation;
    
    fprintf('Starting optimization run %d of %d\n', i, numStarts);
    
    [currentParams, currentFval] = fmincon(@(x)Model_obj(x, p, data), ...
                                            x0_perturbed, ...
                                            [], [], [], [], ... % No linear constraints
                                            LB, UB, ...
                                            [], options);
    
    % Keep the best result from all runs
    if currentFval < bestFval
        bestFval = currentFval;
        bestParams = currentParams;
        fprintf('New best solution found on run %d. Fval = %e\n', i, currentFval);
    end
end

% Assign the best found parameters to be saved
Vmaxs_norm = bestParams;
% Params_norm = bestParams;
fval = bestFval;

save('Vmaxs_norm.mat','Vmaxs_norm');
fprintf('Optimization finished. Best Fval: %e\n', bestFval); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run fitting plot  
run Plot_Fitting.m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print estimated parameters 
disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
disp('Estimated model parameter values');
V_PDH = sprintf('V_PDH = %e', p.pest(p.iPDH)); disp(V_PDH);
V_CITS = sprintf('V_CITS = %e',  p.pest(p.iCITS)); disp(V_CITS);
V_ACON = sprintf('V_ACON = %e',  p.pest(p.iACON)); disp(V_ACON);
V_IDH = sprintf('V_IDH = %e',  p.pest(p.iIDH)); disp(V_IDH);
V_AKGDH = sprintf('V_AKGDH = %e',  p.pest(p.iAKGDH)); disp(V_AKGDH);
V_SCAS = sprintf('V_SCAS = %e',  p.pest(p.iSCAS)); disp(V_SCAS); % 5
V_NDK = sprintf('V_NDK = %e',  p.pest(p.iNDK)); disp(V_NDK);
V_FH = sprintf('V_FH = %e',  p.pest(p.iFH)); disp(V_FH);
V_MDH = sprintf('V_MDH = %e',  p.pest(p.iMDH)); disp(V_MDH);
V_GOT = sprintf('V_GOT = %e',  p.pest(p.iGOT)); disp(V_GOT);
V_CI = sprintf('V_CI = %e',  p.pest(p.iCI)); disp(V_CI); % 10
V_CII_SDH = sprintf('V_CII_SDH = %e',  p.pest(p.iCII)); disp(V_CII_SDH);
V_CIII = sprintf('V_CIII = %e',  p.pest(p.iCIII)); disp(V_CIII);
V_CIV = sprintf('V_CIV = %e',  p.pest(p.iCIV)); disp(V_CIV);
V_CV = sprintf('V_CV = %e',  p.pest(p.iCV)); disp(V_CV); 
V_AK = sprintf('V_AK = %e',  p.pest(p.iAK)); disp(V_AK); 
T_PYR_H = sprintf('T_PYR_H = %e',  p.pest(p.iPYRH)); disp(T_PYR_H); % 15
T_GLU_H = sprintf('T_GLU_H = %e',  p.pest(p.iGLUH)); disp(T_GLU_H);
T_SUC_Pi = sprintf('T_SUC_Pi = %e',  p.pest(p.iDCC1)); disp(T_SUC_Pi);
T_MAL_Pi = sprintf('T_MAL_Pi = %e',  p.pest(p.iDCC2)); disp(T_MAL_Pi);
T_MAL_HCIT = sprintf('T_MAL_HCIT = %e',  p.pest(p.iTCC)); disp(T_MAL_HCIT);
T_MAL_AKG = sprintf('T_MAL_AKG = %e',  p.pest(p.iOME)); disp(T_MAL_AKG); % 20
T_HGHLU_ASP = sprintf('T_HGHLU_ASP = %e',  p.pest(p.iGAE)); disp(T_HGHLU_ASP);
T_ANT = sprintf('T_ANT = %e',  p.pest(p.iANT)); disp(T_ANT);
T_Pi = sprintf('T_PiC = %e',  p.pest(p.iPIC)); disp(T_Pi);
T_HLEAK = sprintf('T_HLEAK = %e',  p.pest(p.iHLEAK)); disp(T_HLEAK); %24
toct=sprintf('time to run the model in min = %d',toc/60); disp(toct); 
 