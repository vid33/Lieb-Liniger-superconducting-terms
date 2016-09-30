function W_out=W_action(t_in, x_in, Q, matF, matr, mass, potential, uu, interaction, D)

        
        R = transpose(reshape(x_in, D, D));
        
        Rkin = Q*R - R*Q;
        
        W_out = -1*( transpose(matF)*R*matr-R*transpose(matF)*matr...
        +(1/(2*mass))*((Rkin*Q - Q*Rkin)*matr+ (Rkin*R - R*Rkin)* matr*R'...
        + (R*R' - R'*R)*Rkin*matr)...
        + potential*R*matr + uu*(matr*R' + R'*matr) ...
        + interaction*(R*R*matr*R' + R'*R*R*matr) );
        
        W_out = reshape(transpose(W_out), D^2, 1);
        
       
end

