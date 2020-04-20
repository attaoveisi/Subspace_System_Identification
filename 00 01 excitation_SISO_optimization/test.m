clc
clear
close all

nFFT = 10;

% Baseband
fmax = 100;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)
df = fmax/nFFT;  % sampling frequency (Hz)
frequency_vector = df:df:fmax;
% frequency_vector = fmax*linspace(0,1,nFFT/2+1);

T_total = nFFT*dt; % duration of the experiment (sec)
time_vector = 0:dt:T_total-dt;
n_inp = 1; % number of input channels

%%
inp_dum = zeros(size(time_vector));
Amp = 1;
ind = 1;
for i = ind:ind:nFFT
    phi_k = -i*(i-1)*pi/nFFT/ind;
    inp_dum = inp_dum + Amp*cos(2*pi*frequency_vector(1,i)*time_vector + phi_k);
end
l_2 = rms(inp_dum);
l_inf = max(abs(inp_dum));

vars = (1:ind:nFFT)*0;
for i = 1:ind:nFFT
    vars(1,i) = -i*(i-1)*pi/nFFT;
end
obj_out = (max(abs(inp_dum)))/l_2;
disp(['initial objective is ', num2str(obj_out)])
disp(['initial crest factor is ', num2str(l_inf/l_2)])
