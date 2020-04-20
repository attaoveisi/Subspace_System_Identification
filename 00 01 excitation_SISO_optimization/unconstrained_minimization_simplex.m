function [x,fval,exitflag,output] = unconstrained_minimization(x0,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimset;
%% Modify options setting
options = optimset(options,'Display', 'off');
options = optimset(options,'MaxFunEvals', MaxFunEvals_Data);
options = optimset(options,'MaxIter', MaxIter_Data);
options = optimset(options,'TolFun', TolFun_Data);
options = optimset(options,'TolX', TolX_Data);
options = optimset(options,'PlotFcn', {  @optimplotfunccount @optimplotfval });
% options = optimset(options,'UseParallel', true);
[x,fval,exitflag,output] = ...
fminsearch(@main_Matlab_objective,x0,options);
