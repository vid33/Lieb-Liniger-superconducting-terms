%%%used by steepest descent/conj. gradient

function [ t_c ] = findOverlapZero_conjGrad_symmetricGauge( R, K, deltaR_orig, deltaK_orig, matr_orig, matl_orig, r_zero, F_zero, Rtangent_orig, D, mass, potential, uu, interaction,...
    t_r, bicg_maxiter, bicg_tolerance, lineSearchTol, maxiter_lineSearch)
    
    t_l = 0;
    
    %preserved by left gauge for tangent vector, so it's the same in all cases.
    matl  = matl_orig; matl_t = transpose(matl);
    
  %  normRtangent_init = real(trace(matl_t*deltaR_orig*matr_orig*deltaR_orig'));
    normRtangent_init = 1;
    
    iterationCounter = 0;

    while 1
        iterationCounter = iterationCounter + 1;

        R_r = R - t_r*deltaR_orig;
        K_r = K - t_r*deltaK_orig;
        K_r = 0.5*(K_r + K_r');

        Rstar_r_l = matl_t\(R_r'*matl_t); 
        
        Q_r = (-1/2)*Rstar_r_l*R_r-1i*(matl_t\K_r);
        Rkin_r = Q_r*R_r - R_r*Q_r;
        [ matr_r, r_zero, final_flag ] = calculate_r_fast( R_r, Q_r, r_zero, D, bicg_maxiter, bicg_tolerance);
        tmp_norm = trace(matl_t*matr_r); matr_r = matr_r/tmp_norm;
        
        [F_r] = calculateFL( R_r, Rkin_r, Q_r, matr_r, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
        F_zero = transpose(F_r);
        matF_r = transpose(reshape(F_r, D, D));
    
        Rstar_r_r = (matr_r*R_r')/matr_r; Rstar_r_l = matl_t\(R_r'*matl_t); 
        Qstar_r_r = (matr_r*Q_r')/matr_r; Qstar_r_l = matl_t\(Q_r'*matl_t); 
        Rtangent_r = matl_t\(transpose(matF_r)*R_r)-R_r*(matl_t\transpose(matF_r))...
            +(1/(2*mass))*(Qstar_r_l*Rkin_r - Rkin_r*Qstar_r_r + R_r*( Rstar_r_l*Rkin_r - Rkin_r*Rstar_r_r)) ...
            +potential*R_r + conj(uu)*(Rstar_r_r + Rstar_r_l)  ...
            + interaction*( Rstar_r_l*(R_r^2) + (R_r^2)*Rstar_r_r);

       % beta = max( real( real(trace(Rtangent_r*matr_r*(Rtangent_r - Rtangent_orig)'))/trace(Rtangent_orig*matr_orig*Rtangent_orig')), 0);
        %fprintf('BETA during linesearch is %d\n', beta);
        
      %  deltaR_r = Rtangent_r + beta*deltaR_orig;        
        
        realOverlap = (1/normRtangent_init)*real(trace(matl_t*Rtangent_r*matr_r*deltaR_orig'));
        
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
            
        
        R_c = R - t_c*deltaR_orig;
        K_c = K - t_c*deltaK_orig;
        K_c = 0.5*(K_c + K_c');

        Rstar_c_l = matl_t\(R_c'*matl_t); 
        
        Q_c = (-1/2)*Rstar_c_l*R_c-1i*(matl_t\K_c);
        
        Rkin_c = Q_c*R_c - R_c*Q_c;
        [ matr_c, r_zero, final_flag ] = calculate_r_fast( R_c, Q_c, r_zero, D, bicg_maxiter, bicg_tolerance);
        tmp_norm = trace(matl_t*matr_c); matr_c = matr_c/tmp_norm;
        
        [F_c] = calculateFL( R_c, Rkin_c, Q_c, matr_c, matl, potential, uu, interaction, F_zero, D, bicg_maxiter, bicg_tolerance );
        F_zero = transpose(F_c);
        matF_c = transpose(reshape(F_c, D, D));

        Rstar_c_r = (matr_c*R_c')/matr_c; Rstar_c_l = matl_t\(R_c'*matl_t); 
        Qstar_c_r = (matr_c*Q_c')/matr_c; Qstar_c_l = matl_t\(Q_c'*matl_t); 
        Rtangent_c = (matl_t\transpose(matF_c)*R_c)-R_c*(matl_t\transpose(matF_c))...
            +(1/(2*mass))*(Qstar_c_l*Rkin_c - Rkin_c*Qstar_c_r + R_c*( Rstar_c_l*Rkin_c - Rkin_c*Rstar_c_r)) ...
            +potential*R_c + conj(uu)*(Rstar_c_r + Rstar_c_l) ...
            + interaction*( Rstar_c_l*(R_c^2) + (R_c^2)*Rstar_c_r);
        
        
       % beta = max( real( real(trace(Rtangent_c*matr_c*(Rtangent_c - Rtangent_orig)'))/trace(Rtangent_orig*matr_orig*Rtangent_orig')), 0);
        %fprintf('BETA is %d\n', beta);
        
       % deltaR_c = Rtangent_c + beta*deltaR_orig;        
        
        realOverlap_c = (1/normRtangent_init)*real(trace(matl_t*Rtangent_c*matr_c*deltaR_orig'));
              
        %fprintf('realOverlap_c is %d\n', realOverlap_c);
        
        if realOverlap_c < 0
            t_r = t_c;
        elseif realOverlap_c > 0
            t_l = t_c;
            
        end
        
    end
    
    
end

