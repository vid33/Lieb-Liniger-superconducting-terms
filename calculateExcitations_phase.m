function [excitation_energies, H, Heigvectors] = calculateExcitations_phase(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, phase, FOR_PLOT)
    
    R1 = phase*R; R2 = R;

    fprintf('Starting calculateExcitations\n');
    fprintf('Preparing necessary matrices...\n');
    excitationPrep_tic = tic;
    
    excitation_energies = zeros(D^2, numel(p)) ;

    Tzero = - kron(Q, eye(D)) - kron( eye(D), conj(Q)) - kron(R, conj(R)) - kron(r,l);
    Tzero_phase = - kron(Q, eye(D)) - kron( eye(D), conj(Q)) - phase*kron(R, conj(R));

    Q1 = Q; Q2 = Q;
    
    R1kin = Q1*R1 - R1*Q1;
    R2kin = Q2*R2 - R2*Q2;
    
    HL = transpose( (1/(2*mass))*(R1kin')*R1kin + potential*(R1')*R1 + uu*R1*R1 + uu*(R1')*(R1') ...
        + interaction*(R1')*(R1')*R1*R1);
    HL = reshape(transpose(HL), 1, D^2);
    HL = HL/Tzero;
    HL = HL - (HL*r)*reshape(eye(D,D),1,D^2); %project on the right
    HL = reshape(HL, D,D);
    

    HR = (1/(2*mass))*R1kin*matr*R1kin' + potential*R1*matr*R1' + uu*R2*R2*matr + uu*matr*(R2')*(R2')...
        + interaction*R1*R1*matr*R1'*R1';
    HR = reshape(transpose(HR), D^2, 1);
    HR = Tzero\HR;
    HR = HR - trace(transpose(reshape(HR, D,D)))*r; %project on the left
    HR = reshape(HR,D,D);
    
    excitationPrep_duration = toc(excitationPrep_tic);
    fprintf('...took %f sec\n', excitationPrep_duration);
    
    if (Tminus_in ~= 0)
        clearvars Tzero;
    end
    
    matrinv_sqrt  = (conj(matr))^(-1/2);
    matr_sqrt  = (conj(matr))^(1/2);

    for k=1: numel(p)
 
        if (FOR_PLOT==1)
            fprintf('%d ', k);
            if (rem(k,10)==0)
                fprintf('\n');
            end
        elseif (FOR_PLOT == 0);
            fprintf('Calculating H...\n'); 
            H_tic = tic;
        end
        
        Q1p = Q1 + R1*R1' + 1i*p(k)*eye(D);
        
        %for large D we save memory by passing Tminus to
        %calculateExcitations. When calculating dispersion relation this is
        %less convenientmatr_sqrt_inv  = (conj(matr))^(-1/2); 
        if (Tminus_in == 0)
            Tminus_phase =  Tzero_phase - 1i*p(k)*eye(D^2);
        else
            Tminus_phase = Tminus_in; 
        end
        
      
        H = kron(R1*(-(1/(2*mass))*Q1p +interaction*eye(D)), matrinv_sqrt*conj(R2)*matr_sqrt) ...
                - (1/(2*mass))*kron(Q1p, matrinv_sqrt*conj(Q2)*matr_sqrt) ...
                +(1/(2*mass))*kron(R1', matrinv_sqrt*conj(Q2)*conj(matr)*transpose(R2)*matrinv_sqrt) ...
            + (1/2)*( kron( ( (1/(2*mass))*(Q1p')*Q1p + potential*eye(D)+interaction*(R1')*R1+HL ), eye(D)) ...
                + kron( ( (1/(2*mass))*R1*(R1')+interaction*eye(D) ), matrinv_sqrt*conj(R2)*conj(matr)*transpose(R2)*matrinv_sqrt) ...
                +kron(eye(D), matrinv_sqrt*( (1/(2*mass))*conj(Q2)*conj(matr)*transpose(Q2) + HR )*matrinv_sqrt ) )...
            +  ((kron( ( (1/(2*mass))*(Q1p')*R1kin + potential*R1 + interaction*(R1')*R1*R1 ...
                + uu*R1' + (HL*R1 - R1*HL)  ), matrinv_sqrt) + uu*kron(eye(D), matrinv_sqrt*conj(R2)) ...
                - kron(  R1*( (1/(2*mass))*R1kin - interaction*R1) ,  matrinv_sqrt*conj(R2)) ...
                - kron( (1/(2*mass))*R1kin, matrinv_sqrt*conj(Q2)))/Tminus_phase)*(kron( eye(D) , conj(R2)*matr_sqrt ) ...
                - kron(R1', matr_sqrt));

        H = (H + H');
        
        if (FOR_PLOT == 0); %end "Calculating H"
            H_duration = toc(H_tic);
            fprintf('...took %f sec\n', H_duration);
        end
        
        fprintf('eig in progress...\n'); 
        eig_tic = tic;
        if (FOR_PLOT == 1)        
            excitation_energies(: , k) = real(eig(H));         
        elseif (FOR_PLOT == 0)
            [Heigvectors, eigvalues] = eig(H);              
            excitation_energies(: , k) = real(diag(eigvalues));
        end
        eig_duration = toc(eig_tic);
        fprintf('... & done after %f sec\n', eig_duration);
    end
      
    if (FOR_PLOT == 1)
        fprintf('\n');
        H = [ ]; Heigvectors =[ ] ;
    end
    
end


        %term4A = kron( transpose(conj(matr)*transpose(R)), eye(D,D)) - kron(hconj(matr), transpose(R));
        %term4A = Tp_plus*term4A;

        %term4A = transpose(term4A);

        %term4A = flip_indices(term4A, D);

        %term4 = term4A*hconj(term3A);
