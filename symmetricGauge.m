
 
    matl_t = transpose(matl);
 [vl,dl]=eig(matl_t);
[dl,indl]=sort(diag(dl),'descend');
vl=vl(:,indl);
    
[vr,dr]=eig(matr);
[dr,indr]=sort(diag(dr),'descend');
vr=vr(:,indr);

Q=diag(1./sqrt(dr))*vr'*Q*vl*diag(1./sqrt(dl));
R=diag(1./sqrt(dr))*vr'*R*vl*diag(1./sqrt(dl));
Z=diag(sqrt(dl))*(vl'*vr)*diag(sqrt(dr));
[UU,SS,VV]=svd(Z);
sqrtS=diag(sqrt(diag(SS)));
Q=sqrtS*VV'*Q*UU*sqrtS;
R=sqrtS*VV'*R*UU*sqrtS;
matr=diag(SS/sqrt(sum(diag(SS).^2)));
%rho(rho<tol/10)=tol/10;
matr=diag(matr);
    
%    break;
%    Q=matr^(-1/2)*Q*matl_t^(-1/2);
%    R=matr^(-1/2)*R*matl_t^(-1/2);
%    Z = matl_t^(1/2)*matr^(1/2);
%    [U,S,V]=svd(Z);
%    Q=S^(1/2)*V'*Q*U*S^(1/2);
%    R= S^(1/2)*V'*R*U*S^(1/2); 

%        matr = S;
    matl = matr; matl_t = transpose(matl);
    
    %if VARY_Q == true
    %           [ matl, matr, l_zero_eigs, r_zero_eigs, l_eigval, r_eigval] = calculate_lr_eigs( R, Q, l_zero_eigs, r_zero_eigs, D, bicg_maxiter, bicg_tolerance);
    %   Q = Q - (r_eigval/2)*eye(D);
    %    [ matl, matr, l_zero_eigs, r_zero_eigs, l_eigval, r_eigval] = calculate_lr_eigs( R, Q, l_zero_eigs, r_zero_eigs, D, bicg_maxiter, bicg_tolerance);
    %   matl_t = transpose(matl);
    %        tmp_norm = trace(matl_t*matr); matr = matr/tmp_norm;

    %end
    

               

    K = 1i*(matl_t*Q + (1/2)*R'*matl_t*R); 
    K = (1/2)*(K + K');
%    K = (1i/2)*(matl_t*Q - Q'*matl_t);

 %   matr_inv = 1./diag(matr); matr_inv = diag(matr_inv); matl_t_inv = matr_inv;

 %   Rstar_l = matl_t_inv*R'*matl_t; 
 %   Q = (-1/2)*Rstar_l*R - 1i*(matl_t_inv*K);
    
    Rkin = Q*R - R*Q;
    
    %    [ matr, r_zero ] = calculate_r_fast( R, Q, r_zero, D, bicg_maxiter, bicg_tolerance);
    %tmp_norm = trace(matl_t*matr); matr = matr/tmp_norm;
    
    l = reshape(transpose(matl), 1, D^2);
    r = reshape(transpose(matr), D^2, 1);