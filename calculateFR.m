function [ FR, flag ] = calculateFR( R, Rkin, Q, matl, matr, potential, uu, interaction, F_zero, D )

    maxiter=20000;
    tolerance=1e-14;

    function out =T_action(in)
        in = transpose(reshape(in, D, D));
        out_tmp=Q*in+in*Q'+R*in*R'+ trace(transpose(matl)*in)*matr;
        out = reshape(transpose(out_tmp), D^2, 1);
    end
    
    rH =  Rkin*matr*Rkin' + potential*R*matr*R' ...
        +uu*(R*R*matr + matr*R'*R') +interaction*R*R*matr*R'*R';

    rH = reshape(transpose(rH), D^2, 1);
    
    [FR ,flag] =bicgstab(@T_action, -rH, tolerance/10, maxiter, [], [], F_zero);
 
    if (flag ~= 0)
        fprintf('NONZERO BICGSTABL FLAG, %d,  IN CALCULATE FR\n', flag);
    end
     
    %FR = FR;

end

