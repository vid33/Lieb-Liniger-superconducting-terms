function [ matr, r_zero_out, final_flag ] = calculate_r_fast( R, Q, r_zero, D, maxiter, tolerance)

    final_flag = 0;
    if nargin < 6
        tolerance=1e-12;
    end
    if nargin < 5
        maxiter=10000;
    end
    
    % max num. of trials before bicgstab decides to quit trying to
    % calculate r to prescribed tolerance.
    errorCount_max = 4;
   
    zero_tmp = zeros(D,D); zero_tmp(1,1) = 1;

    T1 = Q*zero_tmp + zero_tmp*Q' + R*zero_tmp*R';

    T1 = reshape(transpose(T1), D^2, 1);

    T1 = T1(2:end);
    
    [r ,flag, relres, iter] =bicgstab(@Ttilde_action,-T1, tolerance, maxiter, [], [], r_zero(2:end));
    if flag ~= 0
        errorCount = 0; relres_best = relres; r_best = r;
    end
        
    while flag ~= 0
        fprintf('in calculate_r_fast\n');
        fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        fprintf('TRYING bicstab with random r_zero with 4*maxiter!\n');
        
        errorCount = errorCount + 1;
        
        if errorCount >= errorCount_max
            fprintf('Giving up trying to calculate r to prescribed tolerance\n');
            r = r_best; relres = relres_best;
            fprintf('Using best relres: %d\n', relres);
            final_flag = 1;
            break;
        end
      
        r_zero_tmp= cpxrand(D,D)/sqrt(D);
        r_zero_tmp = (1/2)*(r_zero_tmp+r_zero_tmp');
        r_zero_tmp = reshape(transpose(r_zero_tmp), D^2, 1);
        
        [r ,flag, relres, iter] =bicgstab(@Ttilde_action,-T1, tolerance, 4*maxiter, [], [], r_zero_tmp(2:end));
        
        if relres < relres_best
           relres_best = relres; r_best = r; 
        end
        
        fprintf('Is this any better?\n');        
        if flag ~= 0
            fprintf('NO!\n');
            fprintf('At this iteration error in eigenvector_r_fast not corrected\n');
            fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
            %print to file as well
            %fileID = fopen('errorLog', 'w');
            %fprintf(fileID, 'At this iteration error in eigenvector_r_fast not corrected\n');
            %fprintf(fileID, 'Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        else
            fprintf('YES!\n');
            fprintf('Flag: %d, Relres: %d, iter: %d \n', flag, relres, iter);
        end      
    end

    r = [1;r];
 
    r_zero_out = r;
    
    matr = transpose(reshape(r, D, D));
    matr = 0.5*(matr+matr');
    initial_norm = trace(matr);
    matr = (1/initial_norm)*matr;
   
   %r =reshape(transpose(matr), D^2, 1);
   
    function out=Ttilde_action(in)
    
        in = [0; in];
        in = transpose(reshape(in, D,D));
        out=Q*in+in*Q'+R*in*R';
        out = reshape(transpose(out), D^2, 1);
        out = out(2:end);
    end

end