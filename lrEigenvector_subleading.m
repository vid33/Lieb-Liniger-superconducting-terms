function [ l2, r2, matl2, matr2, l2_eigval, r2_eigval] = lrEigenvector_subleading( R, Q, D, eigval_no)

  maxiter=5000;
    tolerance=1e-12;
    
    
    function out=T_action_r(in)
        in = transpose(reshape(in, D, D));
        out=Q*in+in*Q'+R*in*R';
        out = reshape(transpose(out), D^2, 1);
    end

    function out=T_action_l(in)
        in = transpose(reshape(in, D, D));
        out=transpose(Q)*in+in*transpose(Q')+transpose(R)*in*transpose(R');
        out = reshape(transpose(out), D^2, 1);
    end

    opts.disp=0;
    opts.tol=max(tolerance,eps);
    %opts.tol = 1e-3;
    opts.maxit=maxiter;
    opts.isreal=false;
    opts.issym=false;
%    opts.v0=r_zero;

    [r2V,r2_eigvalD, flag]=eigs(@T_action_r,D^2,[],eigval_no,'lr',opts);
    
    r2 = r2V(:,eigval_no);  r2_eigval = r2_eigvalD(eigval_no,eigval_no);
    matr2 = transpose(reshape(r2, D, D));
    
    [l2V,l2_eigvalD,flag]=eigs(@T_action_l,D^2,[],eigval_no,'lr',opts);
    
    l2 = l2V(:,eigval_no);  l2_eigval = l2_eigvalD(eigval_no,eigval_no);
    l2 = transpose(l2);
    matl2 = transpose(reshape(l2, D, D));
    
    norm = l2*r2;
    r2 = r2/norm;
    matr2 = matr2/norm;
    
end