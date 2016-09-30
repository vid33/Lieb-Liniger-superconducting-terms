    [vr,dr]=eig(matr);
    [vl,dl]=eig(transpose(matl));
    Q=dr^(-1/2)*vr'*Q*vl*dl^(-1/2);
    R=dr^(-1/2)*vr'*R*vl*dl^(-1/2);
    Z=dl^(1/2)*(vl'*vr)*dr^(1/2);
    [U,S,V]=svd(Z);
    Q=S^(1)*V'*Q*U;
    R=S^(1)*V'*R*U;
    
    Rkin = Q*R - R*Q;
    
    matr = S/sqrt(sum(diag(S).^2)) ;
    matr_vec = diag(matr); matr = diag(matr_vec).^2;
    matr = matr/trace(matr);
    
    matl = eye(D);
    l = reshape(transpose(matl), 1, D^2);
    r = reshape(transpose(matr), D^2, 1);
    
    K = 1i*(Q + (1/2)*(R')*R);
    K = (1/2)*(K + K');
    Q = (-1/2)*(R')*R -1i*K;
    
    r_zero= cpxrand(D,D)/sqrt(D);
    r_zero = r_zero+r_zero';
    r_zero = reshape(transpose(r_zero), D^2, 1);
    
    [ matr, r_zero ] = calculate_r_fast( R, Q, r_zero, D, 3000, 1e-12);

    matr = matr/trace(matr);

        r = reshape(transpose(matr), D^2, 1);

    
    clearvars S Z U V;