%generate solution files uu -> - uu;

clear;

%dim = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46 48, 50, 52, 54, 56, 58, 60, 64];

D=8;
mass = 0.5;

potential_old =  10;
interaction_old = 0;
uu_old = -5;
Lambda_old = sqrt(20);

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu_old, potential_old, interaction_old);
load(fIn);

dt = 0.000025;

potential =  50;
interaction = 0;
uu = -25;
Lambda = sqrt(20)/sqrt(3);

uu_factor = Lambda/Lambda_old;

R = sqrt(uu_factor)*R;
K = uu_factor*K;
Q = -0.5*(R')*R - 1i*K;
tolerance = 5e-19;
 
Rkin = Q*R - R*Q;

clearvars potential_old interaction_old uu_old Lambda_old uu_factor step;

fOut = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

currentStep=1;
%mainLoop;


