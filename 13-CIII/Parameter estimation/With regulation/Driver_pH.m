% ComplexIII parameter estimation forward reaction pH regulation
% 2CytCo + UQH2 +  2H =  2CytCr + UQ + 4dH
% A - CytCo; B - UQH2; C - CutCr; D -UQ ;

% NOTE: The following code was adopted from the lung model developed
% by Zhang et al 2018 published in Public Library of Science 
% Data used for parameter estimations were obtained from Brandt et al 1997
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024

close all; clear; clc
%% Set initial guesses and bounds for the parameters to be optimized
          
K_H1 = 1.4325200e-07; %Proton binding constant (M)
K_H2 = 1.4325200e-07; %Proton binding constant (M)

p0 = [K_H1,K_H2];  %Binding Constants initial guess

%% Set the optimizer settings and run the optimizer for parameter estimation
% %% Set the optimizer settings and run the optimizer for parameter estimation
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
    lbp = p0/10000; ubp = p0*10000;
    [mpar,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@Error_pH,p0,[],[],[],[],lbp,ubp,[],options);
    
    mpars(i,:) = mpar;
    fvals(i) = fval;
    p0 = temp;
end
[minFval,minIdx] = min(fvals);
mpar = mpars(minIdx,:);
minFval
mpar
save('mpar_pH.txt','mpar','-ascii');

disp('Estimated parameter values'); disp(mpar);
disp('Function at optimization termination'); disp(fval);
disp('Exitflag at optimization termination'); disp(exitflag);
disp('Output at optimization termination'); disp(output);
disp('Lambda at optimization termination'); disp(lambda);
disp('Gradient at optimization termination'); disp(grad);
disp('Hessian at optimization termination'); disp(hessian);

%% Run simulations and plot results
 load mpar_pH.txt
 run Simulations_pH;

