clear;

D=8;
dt=0.002;



%m^2 = 1 lambda = sqrt(20) --> uu=-10/4 potential = 30/4

%gap =0.1 available Ds 4 20 26 32 38

%potential =  -1;
%interaction = 10;
%uu = 1;

mass = 0.5; %m in 1/2m; set this to half.

%interaction = 0;

potential = 0.87141;
interaction = 1;
uu = 1;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d_right.mat', D, D, uu, potential, interaction);

load(fIn);

potential = 0.87140;

K = 0.5*(K + K');
Q = -0.5*(R')*R - 1i*K;

Rkin = Q*R - R*Q;

fOut = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d_right.mat', D, D, uu, potential, interaction);
halt_loop = 0; currentTime = 0;

tic;
tolerance = 5e-17; currentStep =1;
mainLoop_conjGrad;

tolerance = 1e-19; dt = 0.005;
mainLoop_TDVP;
toc;



