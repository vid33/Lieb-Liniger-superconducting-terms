clear;

D = 4;

potential = 10;
interaction = 0;
uu = -5;

fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
load(fIn);

N = 100;
L = 0.05;
position = linspace(0, L, N);

matr_cell = cell(1, N);
matl_cell = cell(1, N);
matFL = cell(1, N);
matFR =  cell(1, N);
testNorm = zeros(1, N);
occupationDensity = zeros(1,N);
eKinDensity = zeros(1,N);
eInteractionDensity = zeros(1,N);
energyDensity = zeros(1, N);
energyDensityR = zeros(1, N);
energyDensityL = zeros(1, N);
entropy = zeros(1, N);
schmidt_sq = zeros(D, N);


dt = 0.0005;

VRsqrt = cpxrand(D,1);
VLsqrt = cpxrand(1,D);

matvr  = kron(VRsqrt', VRsqrt);
matvr = conj(matvr);
vr = reshape(transpose(matvr), D^2, 1);

matvl  = kron(VLsqrt', VLsqrt);
vl = reshape(transpose(matvl), 1, D^2);    

boundaryScript;
