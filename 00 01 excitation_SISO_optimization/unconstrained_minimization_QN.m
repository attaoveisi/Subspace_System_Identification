function [x,fval,exitflag,output,grad,hessian] = unconstrained_minimization_QN(x0,MaxFunctionEvaluations_Data,MaxIterations_Data,OptimalityTolerance_Data,StepTolerance_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('fminunc');
%% Modify options setting
options = optimoptions(options,'Display', 'final-detailed');
options = optimoptions(options,'MaxFunctionEvaluations', MaxFunctionEvaluations_Data);
options = optimoptions(options,'MaxIterations', MaxIterations_Data);
options = optimoptions(options,'OptimalityTolerance', OptimalityTolerance_Data);
options = optimoptions(options,'FunctionTolerance', OptimalityTolerance_Data);
options = optimoptions(options,'StepTolerance', StepTolerance_Data);
options = optimoptions(options,'PlotFcn', {  @optimplotfunccount @optimplotfval });
options = optimoptions(options,'Algorithm', 'quasi-newton');
options = optimoptions(options,'FiniteDifferenceType', 'central');
options = optimoptions(options,'HessUpdate', 'dfp');
options = optimoptions(options,'UseParallel', 'TRUE');
[x,fval,exitflag,output,grad,hessian] = ...
fminunc(@main_Matlab_objective,x0,options);
