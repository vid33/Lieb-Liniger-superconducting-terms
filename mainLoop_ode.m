matl = eye(D);
l = reshape(eye(D), 1, D^2); 

l_zero= cpxrand(D,D)/sqrt(D);
l_zero = l_zero+l_zero';
l_zero = reshape(transpose(l_zero), D^2, 1);

precision = 1e-14;

FIX_GAUGE = true;
DIAGONALIZE_MATR = true;

while 1 
    
    if (currentStep == 1)  
        deltaEnergy = 0;   
        r_zero= cpxrand(D,D)/sqrt(D);
        r_zero = r_zero+r_zero';
        r_zero = reshape(transpose(r_zero), D^2, 1);
        delta_r_zero = r_zero;
        Fzero = cpxrand(D^2, 1);
        delta_Fzero = cpxrand(D^2, 1);
        
        [matr, r_zero_out] = lrEigenvector_fast(R,Q,r_zero, D, precision);
        
        [F, matF ] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, Fzero, D, precision);
        matF_old = matF;

    end
    
    if (currentStep > 1)
              [matr, r_zero_out] = lrEigenvector_fast(R,Q,r_zero, D, precision);
      %  [matr, delta_r_zero] = lrEigenvector_delta(R, Q, matr, delta_r_zero, D, precision);
      r_zero = r_zero_out;
      
    end

    matr = (1/trace(matr))*matr; %r = reshape(matr, 1, D^2);
    
    if (DIAGONALIZE_MATR == true)
   
            [Ur, Sr, Vr] = svd(matr);
            X = inv(Ur);
            matr = X*matr*X';
            R = X*R/X; Q = X*Q/X;  Rkin = Q*R - R*Q;
            K =  (1i/2)*(Q - Q');
            K = 0.5*(K + K');
            matr = real(diag(matr)); matr = diag(matr);
            matr = matr/trace(matr);
            r_zero =reshape(transpose(matr), D^2, 1);
    end
    
        
    calculateEnergy;
    
    if (currentStep > 1)
        deltaEnergy = energyDensity_new - energyDensity_old;
    end 
        
    if (currentStep > 1)  
        [F, matF ] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, Fzero, D, precision);
        F_zero = transpose(F);
        
     %  [ matF, delta_Fzero] = calculateF_delta(R, Rkin, Q, matr, potential, uu, interaction, matF_old, delta_Fzero, D, precision);
     %   matF_old = matF; 
    end
    
%    Rstar = matr*R'/matr;
 %   Rtangent = transpose(matF)*R-R*transpose(matF)...
 %       +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ ((Rkin*R - R*Rkin)*matr*R')/matr...
 %       + (R*R' - R'*R)*Rkin)...
 %       +potential*R + uu*((matr*R')/matr + R')...
 %       +interaction*((R*R*matr*R')/matr+R'*R*R);
%     Rstar = (matr*R')/matr;
%    Rtangent = transpose(matF)*R-R*transpose(matF)...
%        +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ (Rkin*R - R*Rkin)*Rstar...
%        + (R*R' - R'*R)*Rkin)...
%        +potential*R + uu*(Rstar + R')...
%        +interaction*(R*R*Rstar + R'*R*R);
 
    %function M = mass_matrix_ode( t, matr, D)
    
    
    ode_mass_matrix = sparse(transpose(kron(eye(D), matr)));

    %function M = mass_matrix_ode( t, matr, D)
    options_R = odeset('Mass', ode_mass_matrix,  'MassSingular', 'no', 'MStateDependence', 'none', 'MaxStep', dt/4  );
    
    %function W_out=W_action(t_in, x_in, Q, matF, matr, mass, potential, uu, interaction, D)
    [TR, R_vec_new] = ode15s(@(t,x) W_action(t, x, Q, matF, matr, mass, potential, uu, interaction, D), [0 dt/2 dt], reshape(transpose(R), D^2, 1), options_R);
    
%    R_vec_new = deval(R_vec_new_solution, dt);
    
    R_new = transpose(reshape((R_vec_new( end, :)), D, D)); 
   
    Rtangent = (-1/dt)*(R_new - R);
    
    Ktangent = 1i*0.5*(Rtangent'*R-R'*Rtangent);
    Qtangent = -R'*Rtangent;
    
    options_Q = odeset('Mass', ode_mass_matrix,  'MassSingular', 'no', 'MStateDependence', 'none', 'MaxStep', dt/4  );
    [TQ, Q_vec_new] = ode15s(@(t,x) WQ_action(t, x, R, matF, matr, mass, potential, uu, interaction, D), [0 dt/2 dt], reshape(transpose(Q), D^2, 1), options_Q);    
    
%    Q_vec_new = deval(Q_vec_new_solution, dt);
    
    Q_new = transpose(reshape((Q_vec_new(end, :)), D, D)); 
    
    
    normRtangent = trace(Rtangent*matr*Rtangent');
  
    %*<VW | H | Psi(R,Q) > 
    gradOverlap = trace( transpose(matF)*matr*Qtangent' +transpose(matF)*R*matr*Rtangent'...
        +(1/(2*mass))*Rkin*matr*((Rtangent'*Q' - Q'*Rtangent')+(R'*Qtangent' - Qtangent'*R'))...
        +potential*R*matr*Rtangent' ...
        + uu*matr*(R'*Rtangent' + Rtangent'*R') ...
        + interaction*R*R*matr*(R'*Rtangent'+Rtangent'*R') ); 
    
    deltaEnergyLinear = -1*dt*(gradOverlap+conj(gradOverlap));
    
    %    R = R - dt*Rtangent;  
    R = R_new;
    %K = K-dt*Ktangent;
    %K = 0.5*(K + K');
    
    Q = Q_new;
   
    
        %gauge transformation to make matl = eye(D)
    if (FIX_GAUGE == true)

        [ l, matl, l_eigval, l_zero] = eigenvector_l( R, Q, l_zero, D);
        Q = Q -  (1/2)*l_eigval*eye(D);
     
        fprintf('Eigval l is :%d\n', l_eigval);
     
        [Ul, Sl, Vl] = svd(transpose(matl)); 
        X = inv(Ul');
        matl = X'*transpose(matl)*X;  X = inv(X);
        R = (X*R)*inv(X); Q = (X*Q)*inv(X);
        X = matl^(-0.5);
        matl = X'*matl*X;  X = inv(X);
        R = X*R/X; Q = X*Q/X; 
        matl = eye(D); 
    end
    
    Q_tmp = Q;
    K =  (1i/2)*(Q - Q');
    K = 0.5*(K + K');
    
    Q = -0.5*(R')*R-1i*K; 
    Rkin = Q*R - R*Q;
    
    energyDensity_old = energyDensity_new;
    currentStep = currentStep+1;        
    
    if rem(currentStep,1)==0    
        info;
    end
    
    if rem(currentStep,250)==0
        save('data_LLuu/tmp_save.mat');
    end    
        
    if normRtangent< tolerance
        fprintf('Exiting and saving....\n');
        fprintf('Effective energy density is %d, at effective c %d\n',  energyDensity_effective, interaction_effective);
       %save(fOut);
        break;
    end  
     
   % if (currentStep == 250)
   %    toc; break; 
   % end
    
end


