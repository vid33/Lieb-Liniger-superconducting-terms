
%for calculating F and density matrices
bicg_maxiter=10000;
bicg_tolerance=1e-12;

t_c = 0.0001;


%Create a stop button so loop below exits nicely
figh = figure;

global IS_ABORTED;
IS_ABORTED = false;

btn = uicontrol('style', 'pushb', 'string', 'Abort', ...
                'callback', @doAbort);
drawnow;


while 1  
   
    if (currentStep == 1)  
        deltaEnergy = 0;   
         if exist('r_zero', 'var') == 0 
            fprintf('r_zero does not exist. Generating random one. \n');
            r_zero= cpxrand(D,D)/sqrt(D);
            r_zero = r_zero+r_zero';
            r_zero = reshape(transpose(r_zero), D^2, 1);
         end
            
        if exist('F_zero', 'var') == 0 
            fprintf('F_zero does not exist. Generating initial all zeros one. \n');
            F_zero = zeros(D^2, 1);
        end
    end

    if currentStep > 1
        Rtangent_old = Rtangent; matr_old = matr; matl_old = matl; matl_t_old = transpose(matl);
    end
    
  
    calculateEnergy_generalGauge;
    
    if (currentStep > 1)
        deltaEnergy = real(energyDensity_new - energyDensity_old);    
    end                 
    
    Rstar_r = matr*(R'/matr); Rstar_l = (matl_t\R')*matl_t; 
    Qstar_r = matr*(Q'/matr); Qstar_l = (matl_t\Q')*matl_t; 
    Rtangent = (matl_t\transpose(matF))*R-R*(matl_t\transpose(matF))...
        +(1/(2*mass))*(Qstar_l*Rkin - Rkin*Qstar_r + R*( Rstar_l*Rkin - Rkin*Rstar_r))...
        +potential*R + conj(uu)*(Rstar_r + Rstar_l)  ...
        + interaction*(Rstar_l*(R^2) + (R^2)*Rstar_r);
    
    Ktangent = (1i/2)*(Rtangent'*matl_t*R - R'*matl_t*Rtangent);
    Qtangent = -Rstar_l*Rtangent;
    
    normRtangent = trace(matl_t*Rtangent*matr*Rtangent');
  
    gradOverlap = trace( transpose(matF)*matr*Qtangent' +transpose(matF)*R*matr*Rtangent'...
        +(1/(2*mass))*matl_t*Rkin*matr*( Rtangent'*Q' - Q'*Rtangent' + R'*Qtangent' - Qtangent'*R')...
        +potential*matl_t*R*matr*Rtangent' + uu*matl_t*matr*(R'*Rtangent' + Rtangent'*R') ...
        + interaction*matl_t*(R^2)*matr*(R'*Rtangent'+Rtangent'*R') ); 
   
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
         if rem(currentStep, D^2) == 0 || currentStep == 1
            beta = 0;
         else
            beta = max( real(trace(matl_t*Rtangent*matr*(Rtangent - Rtangent_old)'))/real(trace(matl_t_old*Rtangent_old*matr_old*Rtangent_old')), 0);
         end
            
        fprintf('BETA is %d\n', beta);
      
        deltaR = Rtangent + beta*deltaR;
        deltaK = (1i/2)*(deltaR'*matl_t*R-R'*matl_t*deltaR);
        
%            Ktangent = (1i/2)*(Rtangent'*matl_t*R - R'*matl_t*Rtangent);
        
    end
    
    %%%%%BEGIN LINE SEARCH STUFF
    
    t_r_init = 1.1*t_c;
    if currentStep == 1
        t_c = findOverlapZero_generalGauge( R, K, deltaR, deltaK, matr, matl, r_zero, F_zero, Rtangent, D, mass, potential, uu, interaction, t_r_init, bicg_maxiter, bicg_tolerance, 1e-2, 10);    
    elseif currentStep > 1
        t_c = findOverlapZero_conjGrad_generalGauge( R, K, deltaR, deltaK, matr, matl, r_zero, F_zero, Rtangent, D, mass, potential, uu, interaction, t_r_init, bicg_maxiter, bicg_tolerance, 1e-2, 10);
    end
    dt = t_c;
   
   
    deltaEnergyLinear = -1*dt*(gradOverlap+conj(gradOverlap));
    
        
    R = R - dt*deltaR;  
    K = K - dt*deltaK;
    K = 0.5*(K + K');
    
    Rstar_l = matl_t\(R'*matl_t); 
    Q = (-1/2)*Rstar_l*R - 1i*(matl_t\K);
    Rkin = Q*R - R*Q;
    
    %if D ~=1
    [ matr, r_zero, final_flag ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
    tmp_norm = trace(matl_t*matr); matr = matr/tmp_norm;
    %end
    
    if beta == 0 && SYMMETRIC_GAUGE == true
        symmetricGauge;
    end
        
    [F] = calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
    F_zero = transpose(F);
    matF = transpose(reshape(F, D, D));
    
    if final_flag == 1
        fprintf('Right density matrix tolerance not reached. Giving up...\n'); 
        %break;
    end
    
        
    energyDensity_old = energyDensity_new;
    currentStep = currentStep+1;
    
    if rem(currentStep,1)==0      
        info_generalGauge;
    end
    
    if rem(currentStep,250)==0
        save('data_LLuu/tmp_save.mat');
    end  
        
    if real(normRtangent) < tolerance
        [ matr, r_zero, final_flag ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
        tmp_norm = trace(matl_t*matr); matr = matr/tmp_norm;
       
        if SYMMETRIC_GAUGE == true
            symmetricGauge;
        end
            
        [F] = calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
        F_zero = transpose(F);
        matF = transpose(reshape(F, D, D));
        
        fprintf('Exiting and saving....\n');
        fprintf('Effective energy density is %d, at effective c %d\n',  energyDensity_effective, interaction_effective);
        
        close(figh);
       %save(fOut);
        break;
    end  
    
    drawnow;
    if IS_ABORTED == true
        close(figh);
        beta = 0;
        
        [ matr, r_zero, final_flag ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
        tmp_norm = trace(matl_t*matr); matr = matr/tmp_norm;
       
        if SYMMETRIC_GAUGE == true
            symmetricGauge;
        end
            
        [F] = calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
        F_zero = transpose(F);
        matF = transpose(reshape(F, D, D));
        
        break;
    end
    
   
  
end


