% ML estimation of MIMO systems from frequency
% domain data using a subspace method, gradient based search, and the 
% the EM algorithm.
clc
clear
close all

%%
global dsp
dsp = 1;

%%
load('iddata.mat')
ii = 6401;
iddata_act1 = act1_1(1:ii,3)'+act1_1(1:ii,4)'*1i;
iddata_act2 = act2_1(1:ii,3)'+act2_1(1:ii,4)'*1i;
iddata_shaker = shaker_1(1:ii,3)'+shaker_1(1:ii,4)'*1i;
freq_vec = (act1_1(1:ii,2)*2*pi)';

ny = 1;
nu = 1;
nf = ii;
n_x = 12;
fmax = 800;
T = 1/(fmax*2.56);
n = 12;
m = nu;
p = ny;

Y = zeros(ny,nu,nf);
for i = 1:nf
    Y(:,1,i) = iddata_act1(:,i);
%     Y(:,2,i) = iddata_act2(:,i);
%     Y(:,3,i) = iddata_shaker(:,i);
end
%%
z.y = Y; 
w = freq_vec/2/pi;
z.w = w(:); 
z.T = T;

Q = 1*eye(n); R = 1e-1*eye(p); %Initial guess at covariances

mm.A = n; 
mm.op = 'q'; 
mm.T = T; 
mm.type = 'ss';
mm.w = w(:);

oss.alg = 'sid'; % Options for subspace estimation
% oss.lag = round((N-10)/2);
oss.dsp = dsp;

ogn.dsp = dsp; % Options for GN-search estimation
ogn.par = 'ddlc'; 
ogn.cost = 'det'; 
ogn.op = 'q';
ogn.dir = 'trust';
ogn.miter = 100; 
ogn.ngt = 0;

oem.dsp = dsp; % Options for EM-alg estimation
oem.miter = 200; 
oem.optit = 100;
oem.alg = 'em';
oem.stoptol = 1e-6; 

gss = est(z,mm,oss); % Subspace-based estimate
%%
ggn = est(z,gss,ogn); % ML via GN-search estimate starting at sid estimate

%%
% gss.ss.Q = Q; % Reset Q and R matrices to 
% gss.ss.R = R; % Initial values
% gem = est(z,gss,oem); % ML EM search starting at sid estimate
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Display the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp
 data.G = Y; 
 data.w = w(:); 
 data.disp.legend = 'Data';
 
%  showbode(data,gss,gem,ggn);
 showbode(data,gss,ggn);
end

