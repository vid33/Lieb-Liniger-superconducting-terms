matl = eye(D);
l = reshape(eye(D), 1, D^2); 

precision = 1e-14;

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
     %          [matr, r_zero_out] = lrEigenvector_fast(R,Q,r_zero, D, precision);
        [matr, delta_r_zero] = lrEigenvector_delta(R, Q, matr, delta_r_zero, D, precision);
    end

    matr = (1/trace(matr))*matr; %r = reshape(matr, 1, D^2);
        
    calculateEnergy;
    
    if (currentStep > 1)
        deltaEnergy = energyDensity_new - energyDensity_old;
    end 
        
    if (currentStep > 1)  
        [F, matF ] = calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, Fzero, D, precision);
     %  [ matF, delta_Fzero] = calculateF_delta(R, Rkin, Q, matr, potential, uu, interaction, matF_old, delta_Fzero, D, precision);
     %   matF_old = matF; 
    end
    
%    Rstar = matr*R'/matr;
 %   Rtangent = transpose(matF)*R-R*transpose(matF)...
 %       +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ ((Rkin*R - R*Rkin)*matr*R')/matr...
 %       + (R*R' - R'*R)*Rkin)...
 %       +potential*R + uu*((matr*R')/matr + R')...
 %       +interaction*((R*R*matr*R')/matr+R'*R*R);
     Rstar = (matr*R')/matr;
    Rtangent = transpose(matF)*R-R*transpose(matF)...
        +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ (Rkin*R - R*Rkin)*Rstar...
        + (R*R' - R'*R)*Rkin)...
        +potential*R + uu*(Rstar + R')...
        +interaction*(R*R*Rstar + R'*R*R);
    
    solution_r = ode45(@(t,x) T_action_left(t, x, R_spl, Q_spl, D), fliplr(position) , vr, options);

    
    Ktangent = 1i*0.5*(Rtangent'*R-R'*Rtangent);
    Qtangent = -R'*Rtangent;
    
    normRtangent = trace(Rtangent*matr*Rtangent');
  
    %*<VW | H | Psi(R,Q) > 
    gradOverlap = trace( transpose(matF)*matr*Qtangent' +transpose(matF)*R*matr*Rtangent'...
        +(1/(2*mass))*Rkin*matr*((Rtangent'*Q' - Q'*Rtangent')+(R'*Qtangent' - Qtangent'*R'))...
        +potential*R*matr*Rtangent' ...
        + uu*matr*(R'*Rtangent' + Rtangent'*R') ...
        + interaction*R*R*matr*(R'*Rtangent'+Rtangent'*R') ); 
    
    deltaEnergyLinear = -1*dt*(gradOverlap+conj(gradOverlap));
    
    %R_old = R;
    R = R - dt*Rtangent; 
    %R = 0.5*(R+R');
    K = K-dt*Ktangent;
    K = 0.5*(K + K');
    
   % [VR, DR] = eig(R);
   % DR = real(DR);
   % R = VR*DR/VR;
    
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


