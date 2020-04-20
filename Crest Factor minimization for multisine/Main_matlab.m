% Various algorithms were proposed for lowering the multisine crest factor. The following algorithm is slight modification of algorithms proposed in [9] and [10]. It has fast convergence and resulting multisine has good phase randomisation property.
% 
% Algorithm: Low-crest multisine (LCMS) generation.
% 
% Generate uniform random sequence u of length N.
% Transform the sequence to random sign sequence s:
%     for i = 1 : N 
% 
%         if (u(i) > 0.5) 
% 
%             s(i) = 1 
% 
%         else 
% 
%             s(i) = -1 
% 
%     end
% Generate the multisine using the phase values of DFT(s).
% If the number of iterations is less then 10, clip the time series to 85% of the maximum value, else clip to 90% of the maximum value.
% Apply DFT to the clipped time series.
% Generate the new multisine using the phase values obtained in the step 5.
% Calculate the crest factor. If it is higher than objective value go to the step 4.
% Although there is no proof for the convergence of this algorithm, it has never failed in a more than thousand runs with different starting seed factor. The number of iterations grows exponentially when decreasing crest factor. For a multisine sequence of length 32768 it needs only 4-5 iterations to halve the crest factor from 4 to 2, but to get the crest factor < 1.5 it needs several hundreds or more than thousand iterations.
% 
% Phase randomisation properties of this signal can be estimated by testing how it distributes nonlinear distortion in the estimated PIR. The following numerical experiment is done. The input sequence of length 4096 drives 9th order bandpass Butheworth filter (0.1-0.3 fs). The output signal is distorted with 2nd and 3rd order nonlinearities, which are equivalent to 1% of 2nd harmonic distortion and 0.7% of the 3rd harmonic distortion.
% 
% Impulse response is estimated using crosscorrelation and scaled to maximum value 1. Tables I and II show the rms and peak level of the impulse response distortion (in percentage of the maximum value). First 50 values of the impulse response are excluded from distortion calculation to better estimate peak/rms tail distortion. Also, the crest factor of the generated discrete sequence and a crest factor at the filter output are shown.
% 
% These results show that low-crest multisine generates the same level and a distribution of distortions as a true random phase multisine, regardless the crest factor value. MLS sequence generates slightly lower distortion rms level (1-1.5dB), but peak/rms ratio is much higher, especially with 2nd order distortion (10dB higher than with multisine excitation). This high peak/rms ratio disqualifies MLS as a reliable excitation signal for echo detection.
% 
% For all type of an excitation signal the crest factor at the filter output is almost the same. This is a general conclusion for narrow band systems. The benefit of using signals with lower crest factor is significant only for wideband systems, and when it is necessary to use all of the available measurement dynamic range.
% 
% RPMS is the only signal, which do not exhibit a significant crest factor transformation. This fact, constant amplitude and true random phase property fully qualify it as a periodic white noise signal.
% 
% [9] Schoukens, J., Pinelton, R., Ven der Ouderaa, E., and Renneboog, J., "Survey of Excitation Signals for FFT Based Signal Analysers", IEEE Trans. Instrumentation and Measurement, vol. 37, September 1988. [10] Ven der Ouderaa, E., Schoukens, J., and Renneboog, J., "Peak Factor Minimisation of Input and Output Signals of Linear Systems", IEEE Trans. Instrumentation and Measurement, vol. 37, June 1988.

%% script to do a mutlitone low crest factor noise signal
nx = 8192; % fft grid
nf = 200; % number of frequencies
% create a random set of frequencies and gains
ifr = round(nx/2*rand(nf,1))+1;
% gains from -10dB to 0 dB
gains =  10.^(.05*(-10 + 10*rand(nf,1)));

% initialize
phi = 2*pi*rand(nf,1); % random phase
fx0 = zeros(nx/2+1,1);
fx = fx0;
fx(ifr,:) = gains.*exp(1i*phi);
fx0(ifr,:) = gains;

%% main iteration loop
numIter = 20000;
for i = 1:numIter
  % make conjugate symmetric
  fx1 = [fx; conj(fx(end-1:-1:2))];
  % time domain
  x = ifft(fx1);
  % clip at a target crest factor of 2
  xmax = 2*rms(x);
  x( x > xmax)= xmax;
  x( x < -xmax) = -xmax;
  % go back to frequency domain
  fx1 = fft(x);
  % transfer phase
  fx(ifr,:) = gains.*fx1(ifr,:)./abs(fx1(ifr,:));
end
% final time domain signal
fx1 = [fx; conj(fx(end-1:-1:2))];
% time domain
x = ifft(fx1);
fprintf('final crest factor = %f\n', max(abs(x))./rms(x));  