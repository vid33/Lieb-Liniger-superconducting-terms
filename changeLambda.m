clear;

D=32;
dt=0.000001;

potential = 10;
interaction = 0;
uu = -5;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

load(fIn);

Lambda = sqrt(20);

scale_factor = 1.5;

R = sqrt(scale_factor)*R;
K = scale_factor*K;
Lambda = Lambda*scale_factor;

uu = -Lambda^2/4;
potential = Lambda^2/2;

Q = -0.5*(R')*R-1i*K;
Rkin = Q*R - R*Q;

fOut = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

currentStep=1; tolerance = 1e-28;
mainLoop;



