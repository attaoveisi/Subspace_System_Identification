function [x,fval,exitflag,output,population,score] = ga_func(nvars,lb,ub,PopulationSize_Data,CrossoverFraction_Data,MigrationInterval_Data,MigrationFraction_Data,IntCon)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('ga');
%% Modify options setting
options = optimoptions(options,'PopulationSize', PopulationSize_Data);
options = optimoptions(options,'CrossoverFraction', CrossoverFraction_Data);
options = optimoptions(options,'MigrationInterval', MigrationInterval_Data);
options = optimoptions(options,'MigrationFraction', MigrationFraction_Data);
options = optimoptions(options,'CrossoverFcn', @crossoverscattered);
options = optimoptions(options,'MutationFcn', @mutationadaptfeasible);
options = optimoptions(options,'HybridFcn', {  @fmincon [] });
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', {  @gaplotbestf @gaplotbestindiv });
options = optimoptions(options,'UseVectorized', false);
options = optimoptions(options,'UseParallel', true);
[x,fval,exitflag,output,population,score] = ...
ga(@opt_func,nvars,[],[],[],[],lb,ub,[],IntCon,options);
