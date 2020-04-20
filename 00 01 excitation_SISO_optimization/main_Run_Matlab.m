clc
clear
close all
global time_vector
global nFFT
global frequency_vector
global l_2
global ind

%% Fast Fourrier Transformation parameters
N = 400; % 50,100,200,400,800,1600,3200,6400
nFFT = 2.56*N;

% Baseband
fmax = 100;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)
df = fmax/nFFT;  % sampling frequency (Hz)
frequency_vector = df:df:fmax;
% frequency_vector = fmax*linspace(0,1,nFFT/2+1);

T_total = nFFT*dt; % duration of the experiment (sec)
time_vector = 0:dt:T_total-dt;

%% Input signal generation with window and overlap
n_inp = 1; % number of input channels

%%
inp_dum = zeros(size(time_vector));
Amp = 1;
ind = 1;
for i = ind:ind:nFFT
    phi_k = -i*(i-1)*pi/nFFT/ind;
    inp_dum = inp_dum + Amp*cos(2*pi*frequency_vector(1,i)*time_vector + phi_k);
end
l_2 = rms(inp_dum');
l_inf = max(abs(inp_dum));
tic

vars = (1:ind:nFFT)*0;
for i = 1:ind:nFFT
    vars(1,i) = -i*(i-1)*pi/nFFT;
end
obj_out = main_Matlab_objective(vars);
disp(['initial objective is ', num2str(obj_out)])
disp(['initial crest factor is ', num2str(l_inf/l_2)])

%%
% parpool('local')
x0 = vars;
MaxFunEvals_Data = 1e12;
MaxIter_Data = 100000;
TolFun_Data = 1e-5;
TolX_Data = 1e-5;
[x,fval,exitflag,output] = unconstrained_minimization_simplex(x0,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data);
% [x,fval,exitflag,output,grad,hessian] = unconstrained_minimization_QN(x0,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data);
save postotimization_1600

toc
