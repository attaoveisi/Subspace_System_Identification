clc
clear
close all

%% Fast Fourrier Transformation parameters
N = 6400; % 50,100,200,400,800,1600,3200,6400
nFFT = 2.56*N;

% Baseband
fmax = 800;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)
df = fmax/nFFT;  % sampling frequency (Hz)
frequency_vector = df:df:fmax;
% frequency_vector = fmax*linspace(0,1,nFFT/2+1);

T_total = nFFT*dt; % duration of the experiment (sec)
time_vector = 0:dt:T_total-dt;

%% Input signal generation with window and overlap
n_inp = 2; % number of input channels
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
% Creating the input signal (uniformly distributed random signal)
max_inp = 0.5; % specify the maximum amplitude of input
min_inp = -0.5; % specify the maximum amplitude of input
inp_dum = zeros(size(time_vector));
Amp = 1;
for i = 1:nFFT
    phi_k = -i*(i-1)*pi/nFFT;
    inp_dum = inp_dum + Amp*cos(2*pi*frequency_vector(i)*time_vector + phi_k);
end
inp_SISO = inp_dum;
figure(1);
plot(time_vector,inp_SISO)

%%
uniform_win_inp = ones(nFFT,n_inp);
uniform_win_out = ones(nFFT,n_out);

% averaging and overlap properties
n_aver = 1; % number of averaging

% overlap parameters
overlap = 0; 
r_overlap = overlap/100; % rate of overlap (%)
[inp_fft, ~] = pwelch(inp_SISO,uniform_win_out(:,1),overlap,nFFT,df);
inp_fft = inp_fft(2:end)';
inp_fft = [inp_fft fliplr(inp_fft)];
inp_fft_T_dum = zeros(nu,nu,nFFT);
% Hadamard signal creation
for i = 1:nFFT
    inp_fft_T_dum(:,:,i) = inp_fft(i)*T;
end
inp_fft_T = zeros(n_inp,nFFT*nu);
for i = 1:n_inp
    for j = 1:nu
        inp_fft_T(i,(j-1)*nFFT+1:j*nFFT) = squeeze(inp_fft_T_dum(j,i,:));
    end
end
inp_fft_T = inp_fft_T';
inp = zeros(nu*nFFT,n_inp);

for i =1:n_inp
    inp(:,i) = (ifft(inp_fft_T(:,i),'symmetric'));
%     inp(:,i) = inp(:,i)/max(inp(:,i))/2;
end
figure(2)
total_time_vector = 0:dt:((T_total*nu)-dt);
plot(total_time_vector,inp)

%%
T_total = T_total*nu;
T_total = n_aver*T_total-(n_aver-1)*r_overlap*T_total; % total duration of input signal after windowing and averaging (sec)
total_samples = round(n_aver*nFFT-(n_aver-1)*r_overlap*nFFT); % total number of samples
total_time_vector = 0:dt:T_total-dt;
inp = repmat(inp,n_aver,1);
figure(3)
plot(total_time_vector,inp)

%%
load sysc
sample_sys = sysc;
sample_sys = ss(sample_sys.A,sample_sys.B(:,4:5),sample_sys.C,sample_sys.D(:,4:5));

%%
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
        semilogy(freq_after_overlap*1/dt/df,abs(H_2(:,k)))
%         hold on
%         semilogy(f,abs(H_fft(:,k)))
        xlim([0 100])
    end
end

% [sv,w] = sigma(sample_sys);
% figure(4)
% semilogy(w/2/pi,sv(1,:)./sv(2,:))
% xlim([0 100])