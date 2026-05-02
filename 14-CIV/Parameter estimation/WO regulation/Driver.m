
% ComplexIV binding constants' estimation forward reaction w/o regulation
% 2CytCr + 0.5O2 + 2Hm = 2CytCo + H2O + 2dH
% A - CytCr; B - O2 ; C - CytCo
% the flux for CIV reaction in this model was devolped in Audi & Dash lab
% previously(Xiao et al 2018)
% Parameter estimation was done using data from Pannala et al. 2016. and
% Wilson et al 1987 
% Copyright:
% Pardis Taheri, Said Audi, Ranjan Dash 2024

close all; clear; clc

%% Set initial guesses and bounds for the parameters to be optimized

K_A = 7.3375e-04; %CytCr (M)
K_B = 1.5113e-06 ; %O2(M)
Vmaxf = 5.4942e-06; %Maximum forward reaction rate(mmol/min)
 %
   
   
p0 = [K_A, K_B, Vmaxf];  %Binding Constants initial guess
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

