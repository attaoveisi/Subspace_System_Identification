clc
clear
close all
%%
%  PLEASE READ EXERCISE 24 PAGE 42 OF:
%  Johan Schoukens, Rik Pintelon, Yves Rolainauth. Mastering System Identification in 100 Exercises

% In order to run the folder including FDIDENT toolbox must be added.

N = 256;
ExcitationLines = 1 + floor([1:0.1*N]'); % the excited FFT lines
Amp = ones(size(ExcitationLines)); % set the amplitudes
SignalData = fiddata([],Amp,ExcitationLines-1);

UMinCrest = crestmin(SignalData) ; % minimize the crest factor
z = msinprep(UMinCrest,Ndatal) ; % calculates the multisine
uMinCrest = z.input; % put the signal in a vector