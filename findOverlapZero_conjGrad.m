%%%used by steepest descent/conj. gradient

function [ t_c ] = findOverlapZero_conjGrad( R, K, deltaR, deltaK, matr_init, r_zero, F_zero, Rtangent_orig, D, mass, uu, potential, interaction,...
    t_r, bicg_maxiter, bicg_tolerance, lineSearchTol, maxiter_lineSearch)
    
    t_l = 0;
    
%    normRtangent_init = real(trace(deltaR*matr_init*deltaR'));
    normRtangent_init = 1;
    
    iterationCounter = 0;

    while 1
        iterationCounter = iterationCounter + 1;

        R_r = R - t_r*deltaR;
        K_r = K - t_r*deltaK;
        K_r = 0.5*(K_r + K_r');

        K_r = 0.5*(K_r + K_r');
        Q_r = -0.5*(R_r')*R_r-1i*K_r;
        Rkin_r = Q_r*R_r - R_r*Q_r;
        [ matr_r, r_zero ] = calculate_r_fast( R_r, Q_r, r_zero, D, bicg_maxiter, bicg_tolerance);

%        calculateF( R, Rkin, Q, matr, potential, mass, uu, interaction, F_zero, D, maxiter, tolerance )
        
        [F_r] = calculateFL( R_r, Rkin_r, Q_r, matr_r, potential, mass, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
        F_zero = transpose(F_r);
        matF_r = transpose(reshape(F_r, D, D));

        Rstar_r = matr_r*R_r'/matr_r;
        Rtangent_r = transpose(matF_r)*R_r-R_r*transpose(matF_r)...
            +(1/(2*mass))*((Rkin_r*Q_r - Q_r*Rkin_r)+ (Rkin_r*R_r - R_r*Rkin_r)*Rstar_r...
            + (R_r*R_r' - R_r'*R_r)*Rkin_r) + conj(uu)*(Rstar_r + R_r') ...
            +potential*R_r + interaction*(R_r*R_r*Rstar_r + R_r'*R_r*R_r);

       % beta = max( real( real(trace(Rtangent_r*matr_r*(Rtangent_r - Rtangent_orig)'))/trace(Rtangent_orig*matr_orig*Rtangent_orig')), 0);
        %fprintf('BETA during linesearch is %d\n', beta);
        
      %  deltaR_r = Rtangent_r + beta*deltaR_orig;        
        
        realOverlap = (1/normRtangent_init)*real(trace(Rtangent_r*matr_r*deltaR'));
        
        %fprintf('realOverlap is %d\n', realOverlap);
        
        if realOverlap > 0
            t_l = t_r;
            t_r = t_r*2;
        elseif realOverlap < 0
            break;
        end
    
    end 
    iterationCounter = 0;
    while 1
    
        iterationCounter = iterationCounter + 1;
        
        t_c = t_l + (1/2)*(t_r -t_l);
          
        if (t_r - t_l)/t_c < lineSearchTol || iterationCounter > maxiter_lineSearch   
            fprintf('Line search done after %d iterations, dt=%d, realOverlap_c= %d\n', iterationCounter, t_c, realOverlap_c);
            break; 
        end
               
        R_c = R - t_c*deltaR;
        K_c = K - t_c*deltaK;
        K_c = 0.5*(K_c + K_c');

        K_c = 0.5*(K_c + K_c');
        Q_c = -0.5*(R_c')*R_c-1i*K_c;
        Rkin_c = Q_c*R_c - R_c*Q_c;
        [ matr_c, r_zero ] = calculate_r_fast( R_c, Q_c, r_zero, D, bicg_maxiter, bicg_tolerance);

        [F_c] = calculateFL( R_c, Rkin_c, Q_c, matr_c, potential, mass, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
        F_zero = transpose(F_c);
        matF_c = transpose(reshape(F_c, D, D));

        Rstar_c = matr_c*R_c'/matr_c;
        Rtangent_c = transpose(matF_c)*R_c-R_c*transpose(matF_c)...
            +(1/(2*mass))*((Rkin_c*Q_c - Q_c*Rkin_c)+ (Rkin_c*R_c - R_c*Rkin_c)*Rstar_c...
            + (R_c*R_c' - R_c'*R_c)*Rkin_c) + conj(uu)*(Rstar_c + R_c') ...
            +potential*R_c + interaction*(R_c*R_c*Rstar_c + R_c'*R_c*R_c);
        
       % beta = max( real( real(trace(Rtangent_c*matr_c*(Rtangent_c - Rtangent_orig)'))/trace(Rtangent_orig*matr_orig*Rtangent_orig')), 0);
        %fprintf('BETA is %d\n', beta);
        
       % deltaR_c = Rtangent_c + beta*deltaR_orig;        
        
        realOverlap_c = (1/normRtangent_init)*real(trace(Rtangent_c*matr_c*deltaR'));
              
        %fprintf('realOverlap_c is %d\n', realOverlap_c);
        
        if realOverlap_c < 0
            t_r = t_c;
        elseif realOverlap_c > 0
            t_l = t_c;
            
        end
        
    end
    
    
end

