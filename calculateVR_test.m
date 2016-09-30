clear;

potential = 10;
interaction = 0;
uu = -5;

mass = 0.5;

precision = 1e-14;

D=4;
    
fIn = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

load(fIn, 'R', 'Q', 'matr'); matl = eye(D);
Rkin = Q*R - R*Q;
Fzero = cpxrand(D^2, 1);
    
    %function [ F, matF ] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, Fzero, D, precision )

[ FL, matFL] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, Fzero, D, precision);
    %matFL = transpose(reshape(FL, D, D));
[VmatFL, DmatFL] = eig(matFL);
    
[DmatFL, DmatFL_index] = sort(diag(DmatFL));
    
    %for zz=1:D
     
VRsqrt = VmatFL(:, DmatFL_index(1));

    %VRsqrt  = transpose(VRsqrt);
matVR  = kron(VRsqrt', VRsqrt);
matVR = conj(matVR);

VR = reshape(transpose(matVR), D^2, 1);
    
 norm = trace(transpose(matl)*matVR);
    matVR = matVR/norm;
    VR = VR/norm;
    VRsqrt = VRsqrt/sqrt(norm);
   
halt_IDE = 0;
dt = 0.01;  
VRsqrt_test = cpxrand(D,1); counter = 1;
  while halt_IDE == 0
      counter = counter +1;
      
      W_VR = matFL*VRsqrt_test;
      W_VR = W_VR - ((W_VR')*VRsqrt_test)*VRsqrt_test;
        VRsqrt_test = VRsqrt_test - dt*W_VR;
        matVR_test  = kron(VRsqrt_test', VRsqrt_test);
        matVR_test = conj(matVR_test);
        VR_test = reshape(transpose(matVR_test), D^2, 1);
        
        
        norm_tmp = trace(transpose(matl)*matVR_test);
        matVR_test = matVR_test/norm_tmp;
    VR_test = VR_test/norm_tmp;
    VRsqrt_test = VRsqrt_test/sqrt(norm_tmp);
 
      W_VR_norm = trace(conj(kron(W_VR', W_VR)));
      
      fprintf('W_VR_norm is %d\n', W_VR_norm);
      fprintf('l FR overlap is %d\n', trace(transpose(matFL)*matVR_test));
      
      if counter ==100
         break; 
      end
      
  end
  
    
