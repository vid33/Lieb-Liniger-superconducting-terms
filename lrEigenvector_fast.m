function [ matr, r_zero_out ] = lrEigenvector_fast( R, Q, r_zero, D, precision)

    maxiter=10000;

    function out=Ttilde_action(in)
        in = [0; in];
        in = transpose(reshape(in, D,D));
        out=Q*in+in*Q'+R*in*R';
        out = reshape(transpose(out), D^2, 1);
        out = out(2:end);
    end
   
    zero_tmp = zeros(D,D); zero_tmp(1,1) = 1;

    T1 = Q*zero_tmp + zero_tmp*Q' + R*zero_tmp*R';

    T1 = reshape(transpose(T1), D^2, 1);

    T1 = T1(2:end);
    
    [r ,flag] =bicgstab(@Ttilde_action,-T1, max(min(precision/norm(T1), 1e-3)), maxiter, [], [], r_zero(2:end));
    
    if (flag ~= 0)
        fprintf('NONZERO BICGSTABL FLAG, %d,  IN LR EIGENVECTOR\n', flag);
    end

    r = [1;r];
 
    r_zero_out = r;
   
    
   matr = transpose(reshape(r, D, D)); 
   matr = 0.5*(matr+matr');
   

end