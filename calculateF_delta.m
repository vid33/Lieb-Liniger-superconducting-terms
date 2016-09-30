function [ matF,  delta_F ] = calculateF_delta( R, Rkin, Q, matr, potential, uu, interaction, matF_old, delta_Fzero, D, precision )

    maxiter=10000;

    function out =T_action(in)
        in = transpose(reshape(in, D, D));
        out_tmp=transpose(Q)*in+in*conj(Q)+transpose(R)*in*conj(R)+ trace(transpose(matr)*in)*eye(D);
        out = reshape(transpose(out_tmp), D^2, 1);
    end

    
    lH = - (transpose(R)*matF_old*conj(R) + transpose(Q)*matF_old + matF_old*conj(Q)) ...
        -(transpose(Rkin)*conj(Rkin) + potential*transpose(R)*conj(R) ...
        + uu*(transpose(R*R) + conj(R*R)) + interaction*transpose(R*R)*conj(R*R) );
    %- energyDensity_new*eye(D));

    lH = reshape(transpose(lH), D^2, 1);
    
    %fprintf('Norm lH is %d\n', norm(lH));
    
    [delta_F ,flag] =bicgstab(@T_action, lH, max(precision/norm(lH), eps), maxiter, [], [], delta_Fzero);
     
    if (flag ~= 0)
        fprintf('NONZERO BICGSTABL FLAG, %d,  IN CALCULATE F\n', flag);
    end
    
    delta_matF = transpose(reshape(transpose(delta_F), D, D));
    delta_matF = 0.5*(delta_matF + delta_matF');
    delta_F = reshape(transpose(delta_matF), D^2, 1);
    matF = matF_old + delta_matF;
    matF = 0.5*(matF + matF');
    
end

