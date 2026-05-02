% ComplexIV binding constants' estimation forward reaction with pH and mmp
% regulation
% 2CytCr + UQ + 2H -> 2CytCo + UQH2 + 4H(+)

% the flux for CIV reaction in this model was devolped in Audi & Dash lab
% previously(Xiao et al 2018)
% Parameter estimation was done using data from Pannala et al. 2016. 
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024

close all; clear; clc

%% Set initial guesses and bounds for the parameters to be optimized

K_H1 = 7.3375e-04; %CytCr (M)
K_H2 = 1.5113e-06 ; %O2(M)
Vmaxf = 3.5602379e-04; %Maximum forward reaction rate(mmol/min)
 %
   
   
p0 = [K_H1, K_H2, Vmaxf];  %Binding Constants initial guess
%%% Set the optimizer settings and run the optimizer for parameter estimation
options = optimset('fmincon');
options = optimset(options,...
    'TolFun',1e-8,...
    'TolX',1e-8,...
    'Display','iter',...
    'Maxiter',500,...
    'MaxFunEvals',5000);

numIters = 20;
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
    fvals(i) = fval
    p0 = temp;
end

[minFval,minIdx] = min(fvals);
mpar = mpars(minIdx,:);
minFval
mpar

save('mpar.txt','mpar','-ascii');

%disp('Estimated parameter values'); disp(mpar);
%disp('Function at optimization termination'); disp(fval);
% disp('Exitflag at optimization termination'); disp(exitflag);
% disp('Output at optimization termination'); disp(output);
% disp('Lambda at optimization termination'); disp(lambda);
% disp('Gradient at optimization termination'); disp(grad);
% disp('Hessian at optimization termination'); disp(hessian);

%% Run simulations and plot results
 load mpar.txt
 run Simulations;

