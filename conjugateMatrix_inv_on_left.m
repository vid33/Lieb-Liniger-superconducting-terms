function [ Mconj ] = conjugateMatrix_inv_on_left( M, rho )
%CONJUGATEMATRIX_INV_ON_RIGHT  do conjugation rho\M*rho by hand

    rho = diag(rho);
    sizeM = size(M); D = sizeM(1);
    
    Mconj = zeros(D,D);
    
    for kk=1:D
       for mm =1:D 
           if kk == mm
               factor_tmp = 1;
           else
                factor_tmp = rho(mm)/rho(kk);
           end
            Mconj(kk, mm) = M(kk, mm)*factor_tmp;
           
           
       end
    end


end

