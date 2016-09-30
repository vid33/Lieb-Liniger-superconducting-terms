function [ matl, matr, l_zero_eigs_out, r_zero_eigs_out, l_eigval, r_eigval] = calculate_lr_eigs( R, Q, l_zero_eigs, r_zero_eigs, D, maxiter, tolerance)

    function out=T_action_l(in)
        in = (reshape(in, D, D));
        out=transpose(Q)*in+in*conj(Q)+transpose(R)*in*conj(R);
        out = reshape((out), 1, D^2);
    end
    
    
    function out=T_action_r(in)
        in = transpose(reshape(in, D, D));
        out=Q*in+in*Q'+R*in*R';
        out = reshape(transpose(out), D^2, 1);
    end

    opts.disp=0;
    opts.tol=tolerance;
    opts.maxit=maxiter;
    opts.isreal=false;
    opts.issym=false;
    opts.v0=r_zero_eigs;

    [r,r_eigval,flag_r]=eigs(@T_action_r,D^2,[],1,'lr',opts);   
    r_zero_eigs_out = r;
    
    opts.v0=l_zero_eigs;
    [l,l_eigval,flag_l]=eigs(@T_action_l, D^2, [], 1, 'lr', opts);
    l_zero_eigs_out = l;


    if (flag_l ~= 0)
        fprintf('WARNING EIGS MATL FLAG %d AND EIGENVALUE IS %d\n', flag_l, l_eigval);
    end
    
    if (flag_r ~= 0)
        fprintf('WARNING EIGS MATR FLAG %d AND EIGENVALUE IS %d\n', flag_r, r_eigval);
    end
   
    matr = transpose(reshape(r, D, D));
    
    matr = (1/2)*(matr+matr');
   
    matl = transpose(reshape(l, D, D));
    matl = (1/2)*(matl + matl');
    matl = transpose(matl);
   % l = reshape(transpose(matl), 1, D^2);
   
   %r =reshape(transpose(matr), D^2, 1);
   
   %l = reshape(eye(D,D), 1, D^2);

 
    

end