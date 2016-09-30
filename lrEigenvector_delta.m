function [ matr, delta_r] = lrEigenvector_delta( R, Q, matr_old, delta_r_zero, D, precision)

    maxiter=10000;
    
    function out=Ttilde_action(in, R, Q, D, matr_old_tmp)
        in = transpose(reshape(in, D,D));
        out=(Q*in+in*Q'+R*in*R' + trace(in)*matr_old_tmp);
        out = reshape(transpose(out), D^2, 1);
    end
   
    T1 = Q*matr_old + matr_old*Q' + R*matr_old*R';
    T1 = reshape(transpose(T1), D^2, 1);
    
    [delta_r ,flag] =bicgstab(@(in) Ttilde_action(in, R, Q, D, matr_old), -T1, max(min(precision/norm(T1), precision), eps), maxiter, [], [], delta_r_zero);
    
    if (flag ~= 0)
        fprintf('NONZERO BICGSTABL FLAG, %d,  IN LR EIGENVECTOR\n', flag);
    end
    
    delta_matr = transpose(reshape(delta_r, D, D));
    delta_matr = delta_matr - matr_old*trace(delta_matr);
    delta_r = reshape(transpose(delta_r), D^2, 1);
    
    matr = matr_old +delta_matr;    
    matr = 0.5*(matr+matr');
   
end