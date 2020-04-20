clc
clear

%%
nvars = 2;
lb = 17*ones(1,nvars)*3;
ub = 17*ones(1,nvars)*8;
PopulationSize_Data = 6;
CrossoverFraction_Data = 0.95;
MigrationInterval_Data = 2;
MigrationFraction_Data = 0.8;
IntCon = [1 2];
[x,fval,exitflag,output,population,score] = ga_func(nvars,lb,ub,PopulationSize_Data,CrossoverFraction_Data,MigrationInterval_Data,MigrationFraction_Data,IntCon);
