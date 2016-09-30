function [ l, matl, l_eigval, l_zero_out] = eigenvector_l( R, Q, l_zero, D)

    maxiter=10000;
    tolerance=1e-12;

    function out=T_action_l(in)
        in = transpose(reshape(in, D, D));
        out=transpose(Q)*in+in*conj(Q)+transpose(R)*in*conj(R);
        out = reshape(transpose(out), 1, D^2);
    end

    opts.disp=0;
    opts.tol=max(tolerance/100,eps);
    %opts.tol = 1e-3;
    opts.maxit=maxiter;
    opts.isreal=false;
    opts.issym=false;
        opts.v0=l_zero;
        
    [l,l_eigval,flag_l]=eigs(@T_action_l, D^2, [], 1, 'lr', opts);
    
    l_zero_out = l;
    
    %r = 0.5*(r + conj(r));

    if (flag_l ~= 0)
        fprintf('WARNING EIGS FLAG %d AND EIGENVALUE IS %d\n', flag, r_eigval);
    end
   
    matl = transpose(reshape(l, D, D));
    matl = (1/2)*(matl + matl');
    l = reshape(transpose(matl), 1, D^2);
    
    %l_zero_out = transpose(l);
    

end