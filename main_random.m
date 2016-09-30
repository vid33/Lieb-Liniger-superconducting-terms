clear;

D=7;
dt=0.001;

tolerance = 1e-19;

%m^2 = 1 lambda = sqrt(20) --> uu=-10/4 potential = 30/4

%gap =0.1 available Ds 4 20 26 32 38

%potential =  -1;
%interaction = 10;
%uu = 1;

mass = 0.5; %m in 1/2m; set this to half.

%Mass = 10; %unrenormalised mass in \phi^4
%Lambda = sqrt(20); %Cutoff

%potential =  (Mass^2 + Lambda^2)/2;
%uu = (Mass^2 - Lambda^2)/4;
%interaction = 0;

potential = -0.1;
interaction = 5;
uu = 0.1;

R=cpxrand(D,D);

  %[RV, RD] = eig(R);
  %  RD = real(RD); RD  = [1 0 0 0; 0 -1 0 0 ; 0 0 2 0; 0 0 0 -1];
  %  R = RV*RD/RV;

%  R = 0.5*(R+R');
  
K=rand(D,D);
K = 1i*0.5*(K + K');
Q = -0.5*(R')*R - 1i*K;

Rkin = Q*R - R*Q;

matl = eye(D);

fOut = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

currentStep=1; halt_loop = 0; currentTime = 0;

%tic;
mainLoop_TDVP;
tolerance = 1e-25; currentStep =1;
%mainLoop_conjGrad;
%toc;




