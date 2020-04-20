clc
clear
close all

%%
nvars = 2;
lb = 17*ones(1,2);
ub = 120*ones(1,2);
InitialPopulationRange_Data = [30;70];
PopulationSize_Data = 20;
CrossoverFraction_Data = 1.9;
MigrationFraction_Data = 0.9;

[x,fval,exitflag,output,population,score] = ga_opt(nvars,lb,ub,InitialPopulationRange_Data,PopulationSize_Data,CrossoverFraction_Data,MigrationFraction_Data);