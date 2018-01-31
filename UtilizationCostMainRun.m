% News shock script

clear;
clc;
close all;

% Initiallizing parameters
ParameterSetup;

% Steady state solver
SteadyStateSolver;

% Dynare step
dynare UtilizationCostModel.mod noclearall
fprintf('\n')