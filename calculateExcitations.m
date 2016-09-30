function [excitation_energies, H, Heigvectors] = calculateExcitations(p, R, Q, l, r, matr, Tminus_in, D, mass, interaction, potential, uu, FOR_PLOT)

    fprintf('Starting calculateExcitations\n');
    fprintf('Preparing necessary matrices...\n');
    excitationPrep_tic = tic;
    
    excitation_energies = zeros(D^2, numel(p)) ;

    %proj = eye(D^2) - kron(r, l);
    
    Tzero = - kron(Q, eye(D)) - kron( eye(D), conj(Q)) - kron(R, conj(R)) - kron(r,l) ;

    Rkin = Q*R - R*Q;
    
    HL = transpose( (1/(2*mass))*(Rkin')*Rkin + uu*R*R + uu*(R')*(R') ...
        +potential*(R')*R ...
        + interaction*(R')*(R')*R*R);
    HL = reshape(transpose(HL), 1, D^2);
    HL = HL/Tzero;
    HL = HL - (HL*r)*reshape(eye(D,D),1,D^2); %project on the right
    HL = reshape(HL, D,D);
    

    HR = (1/(2*mass))*Rkin*matr*Rkin' + uu*R*R*matr + uu*matr*(R')*(R')...
        + potential*R*matr*R' ...
        + interaction*R*R*matr*R'*R';
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
        
        Qp = Q + R*R' + 1i*p(k)*eye(D,D);
        
        %for large D we save memory by passing Tminus to
        %calculateExcitations. When calculating dispersion relation this is
        %less convenientmatr_sqrt_inv  = (conj(matr))^(-1/2); 
        if (Tminus_in == 0)
            Tminus =  Tzero - 1i*p(k)*eye(D^2, D^2);
        else
            Tminus = Tminus_in; 
        end
        
      
        H = kron(R*(-(1/(2*mass))*Qp +interaction*eye(D,D)), matrinv_sqrt*conj(R)*matr_sqrt) ...
                - (1/(2*mass))*kron(Qp, matrinv_sqrt*conj(Q)*matr_sqrt) ...
                +(1/(2*mass))*kron(R', matrinv_sqrt*conj(Q)*conj(matr)*transpose(R)*matrinv_sqrt) ...
            + (1/2)*( kron( ( (1/(2*mass))*(Qp')*Qp + potential*eye(D,D)+interaction*(R')*R+HL ), eye(D,D)) ...
                + kron( ( (1/(2*mass))*R*(R')+interaction*eye(D,D) ), matrinv_sqrt*conj(R)*conj(matr)*transpose(R)*matrinv_sqrt) ...
                +kron(eye(D), matrinv_sqrt*( (1/(2*mass))*conj(Q)*conj(matr)*transpose(Q) + HR )*matrinv_sqrt ) )...
            +  ( (kron( ( (1/(2*mass))*(Qp')*Rkin + uu*(R') + potential*R + interaction*(R')*R*R ...
                + (HL*R - R*HL)  ), matrinv_sqrt) ...
                + uu*kron(eye(D), matrinv_sqrt*conj(R)) ...
                - kron(  R*( (1/(2*mass))*Rkin - interaction*R) ,  matrinv_sqrt*conj(R)) ...
                - kron( (1/(2*mass))*Rkin, matrinv_sqrt*conj(Q)))/Tminus )*(kron( eye(D) , conj(R)*matr_sqrt ) ...
                - kron(R', matr_sqrt));

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
