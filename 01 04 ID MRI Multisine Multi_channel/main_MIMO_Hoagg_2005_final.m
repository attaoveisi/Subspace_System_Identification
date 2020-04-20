clc
clear
close all

%% Fast Fourrier Transformation parameters
N = 400; % 50,100,200,400,800,1600,3200,6400
nFFT = 2.56*N;

% Baseband
fmax = 100;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)
df = fmax/nFFT;  % sampling frequency (Hz)
frequency_vector = df:df:fmax;

M = ceil(log2(frequency_vector(end)/frequency_vector(1))); % number of frequency segments
dum_Omega = size((2^(M-1))*frequency_vector(1):df:(2^(M))*frequency_vector(1)-df);

Omega = zeros(M,dum_Omega(1,2)); % greatest array needed to assign all Omega segments

N = zeros(M,1); % number of multisines for each frequency segment
for i = 1:M
    omega_beg = (2^(i-1))*frequency_vector(1);
    omega_end = (2^(i))*frequency_vector(1)-df;
    j = size((2^(i-1))*frequency_vector(1):df:(2^(i))*frequency_vector(1)-df,2);
    Omega(i,1:j) = (omega_beg:df:omega_end);
    eta = frequency_vector(end)- (2^(M))*frequency_vector(1); % Page 9 highlighted in red: seems wrong! eta = 0
    eta = df/(1); % the frequency sample of the grid of multisines. select to be half, one fourth,... OR double and...
    N(i) = (omega_end-omega_beg)/eta;
end

%% Input signal generation with window and overlap
n_inp = 3; % number of input channels
n_out = 4; % number of output channels

n_aver = 250; % number of averaging

eta_hat = eta/n_inp; % finer frequency spacing

T_total = 1/eta_hat*M*n_aver; % duration of the experiment (sec)
total_time_vector = 0:dt:T_total-dt;
time_vector = 0:dt:n_aver/eta_hat-dt;
n_samples_i = numel(time_vector);

%% Creating the input signal(sequential multisine)
inp_dum = zeros(size(time_vector));
inp = zeros(n_inp,size(total_time_vector,2));

figure(1);
for k = 1:n_inp
    subplot(n_inp,1,k)
    phase_part = (k-1)*eta_hat;
    for i = 1:M
        for j = 1:N(i,1)
            omega_hat_i = (2^(i-1))*frequency_vector(1);
            Amp = 1/sqrt(N(i,1));
            phi_k = -j*(j-1)*pi/N(i,1);
            inp_dum = inp_dum + Amp*cos(2*pi*(omega_hat_i + eta*(j-1) + phase_part)*time_vector + phi_k);
        end
        inp(k,(i-1)*n_samples_i+1:i*n_samples_i) = inp_dum;
    end
    plot(total_time_vector,inp(k,:))
    pause(0.05)
end

inp = inp';

%%
load sysc
sample_sys = sysc;
if n_inp == 2
    sample_sys = ss(sample_sys.A,sample_sys.B(:,1:2),sample_sys.C,sample_sys.D(:,1:2));
elseif n_inp == 1
    sample_sys = ss(sample_sys.A,sample_sys.B(:,1),sample_sys.C,sample_sys.D(:,1));
else
    load sysc
    sample_sys = ss(sample_sys.A,sample_sys.B(:,1:n_inp),sample_sys.C,sample_sys.D(:,1:n_inp));
end

k = 0;
num_states = sqrt(numel(sample_sys.A));
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        resp(:,k) = squeeze(freqresp(sample_sys(j,i),frequency_vector,'Hz'));
    end 
end
[out_sample,~,~] = lsim(sample_sys,inp,total_time_vector,zeros(num_states,1));

%% Output analysis
out = out_sample;

uniform_win_inp = ones(nFFT,n_inp); % Uniform windowing
uniform_win_out = ones(nFFT,n_out);

overlap = 0; % overlap parameters
r_overlap = overlap/100; % rate of overlap (%)
total_samples = numel(total_time_vector);

if overlap > 0
    n_samples_in_Overlap = floor((overlap*nFFT) / 100); 
    nFrames = n_aver; 
else
    n_samples_in_Overlap= nFFT;
    nFrames=floor(total_samples/n_samples_in_Overlap)-1;
end

S_inp_out = zeros(nFFT/2+1,n_inp*n_out);
k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        [S_inp_out_dum, freq_after_overlap] = cpsd(inp(:,i),out(:,j),uniform_win_inp(:,1),overlap,nFFT,df);
        S_inp_out(:,k) = S_inp_out_dum;
    end
end

S_inp_inp = zeros(nFFT/2+1,n_inp);
for i = 1:n_inp
    [S_inp_inp(:,i), ~] = pwelch(inp(:,i),uniform_win_inp(:,i),overlap,nFFT,df);
end

S_out_out = zeros(nFFT/2+1,n_out);
for j = 1:n_out
    [S_out_out(:,j), ~] = pwelch(out(:,j),uniform_win_out(:,j),overlap,nFFT,df);
end

k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        H_2(:,k) = S_out_out(:,j)./S_inp_out(:,k);
        H_1(:,k) = S_inp_out(:,k)./S_inp_inp(:,i);
    end
end

for i = 1:n_out
    fft_out(:,i) = fft(out(:,i));
end
for i = 1:n_inp
    fft_inp(:,i) = fft(inp(:,i));
end
k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        H_fft(:,k) = fft_out(:,j)./fft_inp(:,i);
    end
end
H_fft = H_fft(1:numel(total_time_vector)/2+1,:);

figure('units','normalized','outerposition',[0 0 1 1])
k = 0;
f = fmax*(0:(numel(total_time_vector)/2))/numel(total_time_vector);
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        subplot(n_inp,n_out,k)
        semilogy(frequency_vector,abs(resp(:,k)))
        hold on
        semilogy(freq_after_overlap*1/dt/df,abs(H_1(:,k)))
        hold on
        semilogy(freq_after_overlap*1/dt/df,abs(H_2(:,k)))
%         hold on
%         semilogy(f,abs(H_fft(:,k)))
        xlim([0 fmax])
    end
end