clc
clear
close all

%% Fast Fourrier Transformation parameters
N = 1600; % 50,100,200,400,800,1600,3200,6400
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
n_inp = 2; % number of input channels
n_out = 2; % number of output channels

%% Hadamard matrix creation for orthogonal input
m = ceil(log2(n_inp)); % 2-79 of Frequency-domain System Identification Book Rik Pintelon, Page 66
nu = 2^m;
H2 = [1 1;1 -1];
if n_inp == 1
    T = 1;
else    
%     Hadamard method
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

    % Dobrowiecki et al. 2006, orthogonal multisine method
%     for pp = 1:n_inp
%         for qq = 1:n_inp
%             T(pp,qq) = (nu^(-1/2))*exp((1i*2*pi*(pp-1)*(qq-1))/nu);
%         end
%     end
end

%%
% averaging and overlap properties
n_aver = 1; % number of averaging

% overlap parameters
overlap = 0; 
r_overlap = overlap/100; % rate of overlap (%)

T_total = n_aver*T_total-(n_aver-1)*r_overlap*T_total; % total duration of input signal after windowing and averaging (sec)
total_samples = round(n_aver*nFFT-(n_aver-1)*r_overlap*nFFT); % total number of samples
total_time_vector = 0:dt:T_total-dt;

% Creating the input signal (uniformly distributed random signal)
max_inp = 1; % specify the maximum amplitude of input
min_inp = -1; % specify the maximum amplitude of input
inp = zeros(total_samples,nu);
for i = 1:total_samples
    phi = -i*(i-1)*pi/total_samples;
    inp(:,:) = inp + cos(2*pi*frequency_vector(1,i)*total_time_vector+phi)'*ones(1,nu);
end
inp_fft = fft(inp);
% Hadamard signal creation
inp_fft = inp_fft*T;
inp = ifft(inp_fft);

% figure(1)
% subplot(2,1,1)
% plot(total_time_vector,inp)
% xlabel('time (sec)')
% ylabel('amplitude')
% subplot(2,1,1)
% n_bins = 500; % number of Histrogram bins
% subplot(2,1,2)
% title('Simulated Histogram of White uniform Noise');
% for i = 1:n_inp
%     [PDF_inp,bins] = hist(inp(:,i),n_bins);
%     hold on
%     bar(bins,PDF_inp/trapz(bins,PDF_inp));
% end
% hold off
% xlabel('bins')
% ylabel('PDF of input signals')
% figure;
% plot(total_time_vector,inp)

%% Test plant
Mass = eye(4);
Stif = [2000 -1000 0 0;-1000 2000 -1000 0;2 -1000 2000 -1000;0 0 -1000 2000];
Damp = blkdiag(20, 5, 0, 0);
if n_inp == 1
    B_inp = [0;2;-1;1];
else
    B_inp = [0,-3;2,0;-1,1;1,0];
end
if n_out == 1
    C_out = [-1,1,-1,1,3,2,-1,2];
else
    C_out = [-1,1,-1,1,5,2,-1,0;1,1,6,-1,0,-2,2,2];
end
% sample_sys = ss([zeros(4) eye(4);-inv(Mass)*Stif -inv(Mass)*Damp],[zeros(4,n_inp);inv(Mass)*B_inp],C_out,zeros(n_out,n_inp));
sample_sys = ss([zeros(4) eye(4);-inv(Mass)*Stif -inv(Mass)*Damp],[zeros(4,n_inp);inv(Mass)*eye(4,n_inp)],C_out,zeros(n_out,n_inp));
% figure(2)
k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
%         subplot(n_inp,n_out,k)
        resp(:,k) = squeeze(freqresp(sample_sys(j,i),frequency_vector,'Hz'));
%         semilogy(frequency_vector,abs(resp))
        %         plot(f,angle(resp)*180/pi);
    end 
end
[out_sample,~,~] = lsim(sample_sys,inp,total_time_vector,zeros(8,1));
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

Win_inp = ones(nFFT,n_inp);
Win_out = ones(nFFT,n_out);

% S_inp_out = zeros(nFFT/2+1,n_inp*n_out);
k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        [S_inp_out_dum, freq_after_overlap] = cpsd(inp(:,i),out(:,j),Win_inp(:,1),overlap,nFFT,df);
        S_inp_out(:,k) = S_inp_out_dum;
    end
end

% S_inp_inp = zeros(nFFT/2+1,n_inp);
for i = 1:n_inp
    [S_inp_inp(:,i), ~] = pwelch(inp(:,i),Win_inp(:,i),overlap,nFFT,df);
end

% S_out_out = zeros(nFFT/2+1,n_out);
for j = 1:n_out
    [S_out_out(:,j), ~] = pwelch(out(:,j),Win_out(:,j),overlap,nFFT,df);
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