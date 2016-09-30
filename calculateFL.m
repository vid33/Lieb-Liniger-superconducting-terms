function [ F] = calculateFL( R, Rkin, Q, matr, matl, potential, uu, interaction, F_zero, D, maxiter, tolerance )

    if nargin < 11
        tolerance=1e-12;
    end
    if nargin < 10
        maxiter=10000;
    end

    % max num. of trials before bicgstab decides to quit trying to
    % calculate F to prescribed tolerance.
    errorCount_max = 4;

    function out =T_action(in)
        in = transpose(reshape(in, D, D));
        out_tmp=transpose(Q)*in+in*conj(Q)+transpose(R)*in*conj(R)+ trace(transpose(matr)*in)*matl;
        out = reshape(transpose(out_tmp), D^2, 1);
    end

    lH = -( transpose(Rkin)*matl*conj(Rkin) + potential*transpose(R)*matl*conj(R) ...
        + uu*transpose(R*R)*matl + conj(uu)*matl*conj(R*R) + interaction*transpose(R*R)*matl*conj(R*R) );
    %lH = -( transpose(Rkin)*conj(Rkin) + potential*transpose(R)*conj(R) ...
    %    +interaction*transpose(R*R)*conj(R*R) );
    %+ energyDensity_new;

    lH = reshape(transpose(lH), D^2, 1);
    
    [F, flag, relres, iter] =bicgstab(@T_action,lH, tolerance, maxiter, [], [], F_zero);
    
    errorCount = 0; relres_best = relres; F_best = F;
    
    while flag ~= 0
        errorCount = errorCount + 1;
        
        if errorCount >= errorCount_max
            fprintf('Giving up trying to calculate F to prescribed tolerance\n');
            F = F_best; relres = relres_best;
            fprintf('Using best relres: %d\n', relres);
            break;
        end
                    
        fprintf('in calculateFL\n');
        fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        fprintf('TRYING bicstab with random F_zero with 4*maxiter!\n');
        
        F_zero_tmp = zeros(D^2, 1);
        
        [F, flag, relres, iter] =bicgstab(@T_action, lH, tolerance, 4*maxiter, [], [], F_zero_tmp);
        
        if relres < relres_best
           relres_best = relres; F_best = F; 
        end
        
        fprintf('Is this any better?\n');        
        if flag ~= 0
            fprintf('NO!\n');
            fprintf('Uncorrected error in calculateLF\n');
            fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
            %print to file as well
            %fileID = fopen('errorLog', 'w');
            %fprintf(fileID, 'Uncorrected error in calculateFL\n');
            %fprintf(fileID, 'Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        else
            fprintf('YES!\n');
            fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        end      
    end
    
    F = transpose(F);

end

