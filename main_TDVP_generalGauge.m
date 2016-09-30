clear;

EXPAND_D = false; %dynamical D expansion; this isn't fully implemented yet.
LOAD_GS = false; %load some ground state solution to use as starting point

VARY_Q =false;  %otherwise vary K..

SYMMETRIC_GAUGE = true;

D=24;
dt=0.0002;
tolerance = 1e-6;

mass = 1/2;
potential = 1;
interaction = 1;
uu = 1;

%for calculating F and density matrices
bicg_maxiter=15000;
bicg_tolerance=1e-13;

if LOAD_GS == true

    fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);
    load(fIn); matl = eye(D);
    
    tolerance = 1e-30;
    
else

    % left gauge
    
    R=cpxrand(D,D);
    K=cpxrand(D,D);
    
    K = (1/2)*(K+ K');
    
    Q = (-1/2)*(R')*R - 1i*(K);
    
    matl = eye(D);
    matl_t = transpose(matl);
    
    r_zero= cpxrand(D,D)/sqrt(D);
    r_zero = r_zero+r_zero';
    r_zero = reshape(transpose(r_zero), D^2, 1);
    F_zero = zeros(D^2, 1);
   
end

[ matr, r_zero ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
tmp_norm = trace(matl_t*matr); matr = matr/tmp_norm;

if SYMMETRIC_GAUGE == true
    symmetricGauge;
end

[F] = calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
F_zero = transpose(F);
matF = transpose(reshape(F, D, D));

if EXPAND_D == true
   
    Qold = Q; Rold = R;
    
    alpha = 10;
    Dold = D; D = D+1;
   
    Rdeformation = zeros(D- Dold);    
    R = [Rold zeros(Dold, D - Dold); zeros(D - Dold, Dold)  Rdeformation ];
    Q = [ Qold zeros(Dold, D - Dold); zeros(D - Dold, Dold) -alpha*eye(D-Dold)];
end


 fOut = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d_generalGauge.mat', D, D, uu, potential, interaction);

if D==1
    l=1; r=1; matl=1; matr=1;
end

currentStep = 1;
currentTime = 0;

%mainLoop;
tic;
mainLoop_TDVP_generalGauge;

%tolerance = 1e-12; SYMMETRIC_GAUGE = false;
%mainLoop_TDVP_generalGauge;

%tolerance = 1e-20; SYMMETRIC_GAUGE = false;

%mainLoop_TDVP_generalGauge;

    
toc;

