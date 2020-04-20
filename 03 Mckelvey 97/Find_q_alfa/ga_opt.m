function [x,fval,exitflag,output,population,score] = ga_opt(nvars,lb,ub,InitialPopulationRange_Data,PopulationSize_Data,CrossoverFraction_Data,MigrationFraction_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('ga');
%% Modify options setting
options = optimoptions(options,'InitialPopulationRange', InitialPopulationRange_Data);
options = optimoptions(options,'PopulationSize', PopulationSize_Data);
options = optimoptions(options,'CrossoverFraction', CrossoverFraction_Data);
options = optimoptions(options,'MigrationFraction', MigrationFraction_Data);
options = optimoptions(options,'FitnessScalingFcn', {  @fitscalingshiftlinear [] });
options = optimoptions(options,'SelectionFcn', @selectionroulette);
options = optimoptions(options,'CrossoverFcn', {  @crossoverheuristic 1.8963 });
options = optimoptions(options,'MutationFcn', @mutationadaptfeasible);
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', {  @gaplotbestf @gaplotbestindiv });
options = optimoptions(options,'UseVectorized', false);
options = optimoptions(options,'UseParallel', true);
[x,fval,exitflag,output,population,score] = ...
ga(@opt_func,nvars,[],[],[],[],lb,ub,[],[],options);
