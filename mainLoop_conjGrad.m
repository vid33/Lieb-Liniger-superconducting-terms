
%for calculating F and density matrices
bicg_maxiter_min=1000;
bicg_tolerance_max=1e-10;

bicg_tolerance_min = 1e-15;
bicg_maxiter_max = 10000;

bicg_maxiter =bicg_maxiter_min;
bicg_tolerance= bicg_tolerance_max;

t_c = 0.0001;

WITNESSED_POS_DELTA_E = 0;

while 1  
   
    if (currentStep == 1)  
        deltaEnergy = 0;   
        r_zero= cpxrand(D,D)/sqrt(D);
        r_zero = r_zero+r_zero';
        r_zero = reshape(transpose(r_zero), D^2, 1);
        F_zero = zeros(D^2, 1);
    end

    if currentStep > 1
        Rtangent_old = Rtangent; matr_old = matr;
    end
    
    if D ~=1
        [ matr, r_zero, final_flag ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
    end
    
    
    if final_flag == 1
        fprintf('Right density matrix tolerance not reached. Giving up...\n'); break;
    end

    
    calculateEnergy;
    
    if (currentStep > 1)
        deltaEnergy = real(energyDensity_new - energyDensity_old);    
    end       
       
        
    [F] = calculateFL( R, Rkin, Q, matr, potential, mass, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
    F_zero = transpose(F);
    matF = transpose(reshape(F, D, D));
    
    Rstar = matr*R'/matr;

    Rtangent = transpose(matF)*R-R*transpose(matF)...
        +(1/(2*mass))*((Rkin*Q - Q*Rkin)+ (Rkin*R - R*Rkin)*Rstar...
        + (R*R' - R'*R)*Rkin)...
        +potential*R + conj(uu)*(Rstar + R')...
        +interaction*(R*R*Rstar + R'*R*R);
    
    normRtangent = trace(Rtangent*matr*Rtangent');
        
    Ktangent = 1i*0.5*(Rtangent'*R-R'*Rtangent);
    Qtangent = -R'*Rtangent;
    
  
    gradOverlap = trace( transpose(matF)*matr*Qtangent' +transpose(matF)*R*matr*Rtangent'...
        +(1/(2*mass))*Rkin*matr*((Rtangent'*Q' - Q'*Rtangent')+(R'*Qtangent' - Qtangent'*R'))...
        +potential*R*matr*Rtangent' ...
        + conj(uu)*matr*(R'*Rtangent' + Rtangent'*R') ...
        + interaction*R*R*matr*(R'*Rtangent'+Rtangent'*R') ); 
   
    if currentStep == 1
        deltaR = Rtangent;
        deltaK = Ktangent;
        deltaQ = Qtangent;
        Rtangent_old = Rtangent;
        
    elseif currentStep > 1
        %fletcher-reeves
   %     beta = real(trace(Rtangent*matr*Rtangent')/trace(Rtangent_old*matr_old*Rtangent_old'));
        
         %polak-ribiere
         
         %restart conjGrad every D^2 steps
         if rem(currentStep, D^2) == 0;
             beta = 0;
         else
            beta = max( real(trace(Rtangent*matr*(Rtangent - Rtangent_old)'))/real(trace(Rtangent_old*matr_old*Rtangent_old')), 0);
         end
            
        fprintf('BETA is %d\n', beta);
      
        deltaR = Rtangent + beta*deltaR;
        deltaK = 1i*0.5*(deltaR'*R-R'*deltaR);
    end
    
    %%%%%BEGIN LINE SEARCH STUFF
    
    t_r_init = 1.1*t_c;
    if currentStep == 1
        t_c = findOverlapZero( R, K, deltaR, deltaK, matr, r_zero, F_zero, Rtangent, D, mass, uu, potential, interaction, t_r_init, bicg_maxiter, bicg_tolerance, 1e-2, 10);    
    elseif currentStep > 1
        t_c = findOverlapZero_conjGrad( R, K, deltaR, deltaK, matr, r_zero, F_zero, Rtangent, D, mass, uu, potential, interaction, t_r_init, bicg_maxiter, bicg_tolerance, 1e-2, 10);

    end
    dt = t_c;
   
   
    deltaEnergyLinear = -1*dt*(gradOverlap+conj(gradOverlap));
    
        
    R = R - dt*deltaR;  
    K = K - dt*deltaK;
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
        
    if real(normRtangent) < tolerance
        fprintf('Exiting and saving....\n');
        fprintf('Effective energy density is %d, at effective c %d\n',  energyDensity_effective, interaction_effective);
       %save(fOut);
        break;
    end  
   
  
end


