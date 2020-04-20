clc
clear
close all

%% Fast Fourrier Transformation parameters
N = 6400; % 50,100,200,400,800,1600,3200,6400
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
n_out = 6; % number of output channels

%% Hadamard matrix creation for orthogonal input
m = ceil(log2(n_inp)); % 2-79 of Frequency-domain System Identification Book Rik Pintelon, Page 66
if m == 0
    m = m+1;
end
nu = 2^m;
H2 = [1 1;1 -1];
if m == 1
    H2_dum = H2;
else
    for i = 1:m-1
        if i == 1
            H2_dum = H2;
        end
        H2_dum = kron(H2,H2_dum); 
    end
end
H2m = H2_dum; 
T = 1/sqrt(nu)*H2m;

%%
% Windowing (Hanning window): w(n)=0.5(1?cos(2?n/N)), 0?n?N L=N+1
% w = hann(L,'sflag') returns an L-point Hann window using the window sampling specified by 'sflag',
% which can be either 'periodic' or 'symmetric' (the default). 
% The 'periodic' flag is useful for DFT/FFT purposes, such as in spectral analysis.

hann_win_inp = zeros(nFFT,n_inp);
for i = 1:n_inp
    hann_win_inp(:,i) = hann(nFFT,'periodic');
end

hann_win_out = zeros(nFFT,n_inp);
for i = 1:n_out
    hann_win_out(:,i) = hann(nFFT,'periodic');
end

% averaging and overlap properties
n_aver = 100; % number of averaging

% overlap parameters
overlap = 75; 
r_overlap = overlap/100; % rate of overlap (%)

T_total = n_aver*T_total-(n_aver-1)*r_overlap*T_total; % total duration of input signal after windowing and averaging (sec)
total_samples = round(n_aver*nFFT-(n_aver-1)*r_overlap*nFFT); % total number of samples
total_time_vector = 0:dt:T_total-dt;

% Creating the input signal (uniformly distributed random signal)
max_inp = 0.5; % specify the maximum amplitude of input
min_inp = -0.5; % specify the maximum amplitude of input
% inp = max_inp*ones(size(randn(total_samples,n_inp))) + (max_inp-min_inp)*randn(total_samples,n_inp);

% 'rgs' — Gives a random, Gaussian signal.
% 'rbs' — Gives a random, binary signal. This is the default.
% 'prbs' — Gives a pseudorandom, binary signal.
inp_SISO = idinput([numel(total_time_vector) 1],'rgs',[0 1],[min_inp max_inp]);

%%
inp_fft = fft(inp_SISO);
inp_fft_T_dum = zeros(nu,nu,numel(total_time_vector));
% Hadamard signal creation
for i = 1:numel(total_time_vector)
    inp_fft_T_dum(:,:,i) = inp_fft(i)*T;
end
inp_fft_T = zeros(n_inp,numel(total_time_vector)*nu);
for i = 1:n_inp
    for j = 1:nu
        inp_fft_T(i,(j-1)*numel(total_time_vector)+1:j*numel(total_time_vector)) = squeeze(inp_fft_T_dum(j,i,:));
    end
end
inp_fft_T = inp_fft_T';
inp = zeros(nu*numel(total_time_vector),n_inp);
for i =1:n_inp
    if i == 1
        inp(:,i) = (ifft(inp_fft_T(:,i),'symmetric'));
    else
        inp(:,i) = (ifft(inp_fft_T(:,i),'symmetric'));
    end
end
%%
% revise the time vector
total_time_vector = 0:dt:((T_total*nu)-dt);
figure(1)
subplot(2,1,1)
plot(total_time_vector,inp)
xlabel('time (sec)')
ylabel('amplitude')
subplot(2,1,1)
n_bins = 500; % number of Histrogram bins
subplot(2,1,2)
title('Simulated Histogram of White uniform Noise');
for i = 1:n_inp
    [PDF_inp,bins] = hist(inp(:,i),n_bins);
    hold on
    bar(bins,PDF_inp/trapz(bins,PDF_inp));
end
hold off
xlabel('bins')
ylabel('PDF of input signals')
figure;
plot(total_time_vector,inp)

%% Test plant
load sysc
sample_sys = sysc;
% sample_sys = ss(sample_sys.A,sample_sys.B(:,4:5),sample_sys.C,sample_sys.D(:,4:5));
sample_sys = ss(sample_sys.A,sample_sys.B(:,4),sample_sys.C,sample_sys.D(:,4));

% figure(2)
k = 0;
num_states = sqrt(numel(sample_sys.A));
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
%         subplot(n_inp,n_out,k)
        resp(:,k) = squeeze(freqresp(sample_sys(j,i),frequency_vector,'Hz'));
%         semilogy(frequency_vector,abs(resp))
        %         plot(f,angle(resp)*180/pi);
    end 
end
[out_sample,~,~] = lsim(sample_sys,inp,total_time_vector,zeros(num_states,1));
% figure(3)
% plot(total_time_vector,out_sample)

%% Output analysis
out = out_sample;

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
        [S_inp_out_dum, freq_after_overlap] = cpsd(inp(:,i),out(:,j),hann_win_inp(:,1),overlap,nFFT,df);
        S_inp_out(:,k) = S_inp_out_dum;
    end
end

S_inp_inp = zeros(nFFT/2+1,n_inp);
for i = 1:n_inp
    [S_inp_inp(:,i), ~] = pwelch(inp(:,i),hann_win_inp(:,i),overlap,nFFT,df);
end

S_out_out = zeros(nFFT/2+1,n_out);
for j = 1:n_out
    [S_out_out(:,j), ~] = pwelch(out(:,j),hann_win_out(:,j),overlap,nFFT,df);
end

k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        H_2(:,k) = S_out_out(:,j)./S_inp_out(:,k);
    end
end

figure(3)
k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        subplot(n_inp,n_out,k)
        semilogy(frequency_vector,abs(resp(:,k)))
        hold on
        semilogy(freq_after_overlap*1/dt/df,abs(H_2(:,k)))
        xlim([0 100])
    end
end

% [sv,w] = sigma(sample_sys);
% figure(4)
% semilogy(w/2/pi,sv(1,:)./sv(2,:))
% xlim([0 100])