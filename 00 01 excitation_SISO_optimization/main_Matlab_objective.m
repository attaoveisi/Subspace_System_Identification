function obj_out = main_Matlab_objective(vars)

global time_vector
global nFFT
global frequency_vector
global l_2
global ind

%%%  Very important %%%
p = 8;

inp_dum = zeros(size(time_vector));
Amp = 1;
phi_k = vars;
for i = 1:ind:nFFT
    inp_dum = inp_dum + Amp*cos(2*pi*frequency_vector(i)*time_vector + phi_k(i));
end

% inp_dum = zeros(nFFT,size(time_vector,2));
% Amp = 1;
% phi_k = vars;
% for i = 1:nFFT
%     inp_dum(i,:) = Amp*cos(2*pi*frequency_vector(i)*time_vector + phi_k(i));
% end
% inp_dum = sum(inp_dum);

% inp_dum = zeros(nFFT,size(time_vector,2));
% Amp = 1;
% phi_k = vars;
% for i = 1:nFFT
%     inp_dum(i,:) =  2*pi*frequency_vector(i)*time_vector + phi_k(i);
% end
% inp_dum = sum(Amp*cos(inp_dum));

% obj_out = (1/numel(time_vector)*sum(abs(inp_dum).^p))/l_2;
obj_out = (max(abs(inp_dum)))/l_2;