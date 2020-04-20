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
% frequency_vector = fmax*linspace(0,1,nFFT/2+1);
ww = frequency_vector*2*pi;
TT = 5/max(ww);
wk = 2*atan(TT*ww/2);

T_total = nFFT*dt; % duration of the experiment (sec)
time_vector = 0:dt:T_total-dt;

%% Input signal generation with window and overlap
n_inp = 7; % number of input channels
n_out = 6; % number of output channels

%% Dobrowiecki matrix creation for orthogonal input
T = zeros(n_inp,n_inp);
for i = 1:n_inp
    for j = 1:n_inp
        T(i,j) = (n_inp^(-0.5))*exp(1i*2*pi*(i-1)*(j-1)/n_inp);
    end
end

%%
% Creating the input signal (uniformly distributed random signal)
max_inp = 0.5; % specify the maximum amplitude of input
min_inp = -0.5; % specify the maximum amplitude of input
inp_dum = zeros(size(time_vector));
num_sum = fmax;
stepp = fmax/N;
Amp = 1/sqrt(numel(1:stepp:num_sum));
frequency_vector1 = 1:stepp:fmax;
ii = 1;
for i = 1:stepp:num_sum
    phi_k = -ii*(ii-1)*pi/num_sum;
%     phi_k = x(i);
    inp_dum = inp_dum + Amp*cos(2*pi*frequency_vector1(ii)*time_vector + phi_k);
    ii = ii+1;
end
inp_SISO = inp_dum;
figure(1);
subplot(2,1,1)
plot(time_vector,inp_SISO)
subplot(2,1,2)
fff = fs*(0:(nFFT/2))/nFFT;
% inp_SISO(1,1025:end) = inp_SISO(1,1:1024);
inp_FFT = fft(inp_SISO);
plot(fff,20*log10(abs(inp_FFT(1,1:nFFT/2+1))))

%%
inp_fft = fft(inp_SISO);
inp_fft_T_dum = zeros(n_inp,n_inp,nFFT);
% Hadamard signal creation
for i = 1:nFFT
    inp_fft_T_dum(:,:,i) = inp_fft(i)*T;
end
inp = zeros(n_inp*nFFT,n_inp);
for j = 1:n_inp
    for i = 1:n_inp
%         inp(i,(j-1)*nFFT+1:j*nFFT) = ifft(squeeze(inp_fft_T_dum(i,j,:)))';
        inp((j-1)*nFFT+1:j*nFFT,i) = ifft(squeeze(inp_fft_T_dum(i,j,:)),'symmetric');
    end
end

figure(2)
subplot(2,1,1)
total_time_vector = 0:dt:((T_total*n_inp)-dt);
plot(total_time_vector,inp)
subplot(2,1,2)
fff = fs*(0:(nFFT/2))/nFFT;
hold on
for i = 1:n_inp
    inp_FFT = abs(fft((inp(:,i))'));
    hold on
    plot(fff,20*log10((inp_FFT(1,1:nFFT/2+1))),'*')
end

%%
uniform_win_inp = ones(nFFT,n_inp);
uniform_win_out = ones(nFFT,n_out);

% averaging and overlap properties
n_aver = 1; % number of averaging

% overlap parameters
overlap = 0; 
r_overlap = overlap/100; % rate of overlap (%)
T_total = T_total*n_inp;
T_total = n_aver*T_total-(n_aver-1)*r_overlap*T_total; % total duration of input signal after windowing and averaging (sec)
total_samples = round(n_aver*nFFT-(n_aver-1)*r_overlap*nFFT); % total number of samples
total_time_vector = 0:dt:T_total-dt;
inp = repmat(inp,n_aver,1);
% figure(3)
% plot(total_time_vector,inp)

%%
load sysc
sample_sys = sysc;
% sample_sys = ss(sample_sys.A,sample_sys.B(:,1:2),sample_sys.C,sample_sys.D(:,1:2));
%%
% i = 1;
% while i
%     sysc = rss(14,6,7);
%     sample_sys  = sysc;
% %     sys.A = rand(14,14) + [-0.312459745646968,58.6755716242942,-3.32970939842960,10.3001198757494,-6.95387719278478,-3.36814766828286,-1.21486316836128,-3.56545637290002,1.75564668058367,-1.30668435176006,2.69923478036896,0.994098389511011,-0.339309909093562,0.131241118634470;-56.9226724070856,-2.75113184971002,6.66402072060328,0.340905151945994,3.38917171329582,23.2194873322224,1.83037894267740,3.67995347699375,3.06552440400818,-2.21056858336509,-6.39178071982871,-4.26800925655225,1.25178756966999,0.0856989227708177;2.27059054207520,-2.86378857809285,-2.39928036339358,-91.0922217774496,-50.2647625309475,-45.2274829495707,-0.480167334049433,0.438103167323734,6.94528951098814,8.02855997653458,1.50196034098208,1.17353369457493,-0.0604880356552731,0.251559976404786;-8.52202820906988,-6.10381802878106,92.3953302206382,-4.65829173669476,-6.62790878547604,192.592133159742,-2.47132108526682,13.1281484997033,2.33061854546082,-2.02226918844728,-9.27795934270479,-4.60935560511868,1.43064371741356,-0.196358628355626;5.43107685534713,1.27706813347834,45.9816043439077,11.4919101376484,-1.92453700185895,-150.813064205593,0.259015713237698,-3.36249298155650,-5.25790555794321,-1.78151950310267,5.38951136906085,3.49744450286000,-1.69086108534427,0.180627921804192;1.20916351112345,-15.2957511580614,41.1608984759840,-182.159977647126,144.151510485498,-5.89816234851918,-30.7906087108525,1.26261206413858,-3.96048614873784,5.76698987677539,5.24809114362927,8.99399522884542,-3.05230259010260,0.208494668207586;0.641520180823860,-0.581774347381974,-0.752406772651269,4.70617515807350,-1.66121993924186,28.3710736029170,-2.93573445227263,-163.786248215560,1.64422730203414,14.4874903585412,-3.06973058933075,-8.83632290707219,0.124637759337413,0.960109991455750;2.00095537741955,0.567136525731416,-3.30045542779170,-7.38997677973473,-0.561818311601700,-8.47029695846791,158.027273754072,-3.54434194091480,-24.0167729101540,8.87157905904816,-2.62820867070228,2.13587947920612,-2.81920390874451,0.826397190946874;-2.28717831671985,-3.43627237240055,-7.23515519385282,-1.45514934033785,3.96097211339719,2.76000493920257,-6.51842734214816,19.1007507544078,-7.70823576790647,-1.85165475864283,8.49277415624101,-1.15630189083034,-5.53070327894554,2.77308184992660;1.62256847137151,2.06761940740116,-6.39766044203789,-0.820106788227726,2.75811641287774,-3.18021260091163,-11.3408874378579,-6.78442848680085,1.69331543256094,-9.35588873243869,-1.62444277326928,-0.627845712389132,-0.121163482663899,0.143757560683648;-0.453345513880598,-1.40783541828260,2.04180660869143,-0.668786133169398,1.87239810362728,6.09328210889553,-1.44286743426222,7.04364101258553,-4.31222855157339,-0.392114125133970,-22.3742028433171,-186.207075689043,46.2207300229101,-0.714320122730896;-0.0382791428405690,1.84539593327961,-1.83120708894555,1.03130392765232,-1.77745160397810,-7.38547092906018,3.01968108484277,-5.09807064224737,0.797666266425589,0.661070245739320,177.719240523915,-11.2684272915322,2.27396484281817,0.407911224374075;-0.170232319350042,-0.707506781620471,0.0362610052826425,-0.252801528000997,0.590764462355496,2.04683816313453,-1.31050150889338,1.53494269397176,-1.48272899696038,-0.316481070731215,-44.5439853251676,4.00458064625061,-2.35932720087457,6.42912196800845;-0.0716260780736722,-0.102669791570725,-0.174773033429599,0.0282082849242810,0.0252887640646545,0.00961633119060952,-0.344954823207545,-0.00172763288042097,-1.01515748268251,-0.0570264517230666,-0.670202946476633,1.06841075856952,-5.63472580522611,-2.35825242790048];
%     if rank(ctrb(sysc)) == 14
%         if rank(obsv(sysc)) == 14
%             i = 0;
%         end
%     end
% end
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
%         hold on
%         semilogy(f,abs(H_fft(:,k)))
        xlim([10 fmax])
    end
end

% [sv,w] = sigma(sample_sys);
% figure(4)
% semilogy(w/2/pi,sv(1,:)./sv(2,:))
% xlim([0 100])