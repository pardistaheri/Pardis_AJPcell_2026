% Citrate synthase parameter estimation with ATP, ADP, AMP, and SCOA regulation
% ACOA + OXA = COA + CIT + 2H+
% A - ACOA; B - OXA; C - COA; D - CIT;
%
% NOTE: The following regulation was adopted from the citrate synthase model developed
% by Beard et al 2008 published in Public Library of Science 
% Data used for parameter estimations were obtained from Shepherd et al
% 1969, and Smith et al 1971
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024
% Medical College of Wisconsin and Marquette University joint program in Biomedical Engineering

clear all;
close all;
clc;

% Set initial guesses and bounds for the parameters to be optimized
K_ATP = 37.3*1e-5; %ATP binding constant in isolated CITS(M)
K_ADP = 135.4*1e-4; %ADP binding constant in isolated CITS(M)
K_AMP = 992.3*1e-4; %AMP binding constant in isolated CITS(M)
K_SCOA = 74.1e-06; %SCoA binding constant (M)


p0 = [K_ATP, K_ADP, K_AMP,K_SCOA];  %Binding Constants initial guess
%%% Set the optimizer settings and run the optimizer for parameter estimation
options = optimset('fmincon');
options = optimset(options,...
    'TolFun',1e-8,...
    'TolX',1e-8,...
    'Display','iter',...
    'Maxiter',500,...
    'MaxFunEvals',5000);

numIters = 100;
mpars = []; fvals = [];
randLb = 0.5; randUb = 1.5;
numPars = length(p0);

for i = 1:numIters
    temp = p0;
    randNum = (randLb + (randUb-randLb)*rand(1,numPars));
    %Generates random number b/w 0.5 and 1.5 effectively allowing p0 to
    %vary between 0.5p0 and 1.5p0 (a maximal increase or decrease of 50%) 

    p0 = p0.*randNum;
    lbp = p0/100; ubp = p0*100;
    [mpar,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@Error,p0,[],[],[],[],lbp,ubp,[],options);
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
    p0 = temp;
end

[minFval,minIdx] = min(fvals);
mpar = mpars(minIdx,:);
minFval
mpar

save('mpar.txt','mpar','-ascii');

% disp('Estimated parameter values'); disp(mpar);
% disp('Function at optimization termination'); disp(fval);
% disp('Exitflag at optimization termination'); disp(exitflag);
% disp('Output at optimization termination'); disp(output);
% disp('Lambda at optimization termination'); disp(lambda);
% disp('Gradient at optimization termination'); disp(grad);
% disp('Hessian at optimization termination'); disp(hessian);

%% Run simulations and plot results
 load mpar.txt
 run Simulations

