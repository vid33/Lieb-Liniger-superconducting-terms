function out=F_action_left_const(t_in, F_in, R, Q, Rkin, mass, potential, uu, interaction, D, solution_r)
        F_in = transpose(reshape(F_in, D, D));
       % F_in = conj(F_in);
        r_tmp = deval(solution_r, t_in);
          
        matr_tmp = transpose(reshape((r_tmp), D, D));
        %matr_tmp = conj(matr_tmp);
        out_tmp= Q*F_in  + F_in*Q'  + R*F_in*R' ...
                + (1/(2*mass))*Rkin*matr_tmp*Rkin' + uu*(R*R*matr_tmp  + matr_tmp*(R')*(R')) ...
                + potential*R*matr_tmp*R' + interaction*R*R*matr_tmp*(R')*(R');
            %- energyDensity*matr_tmp;
        out = -reshape(transpose(out_tmp), D^2, 1);

end

