function out=F_action_right_const(t_in, F_in, R, Q, Rkin, mass, potential, uu, interaction, D, solution_l)
        F_in = transpose(reshape(F_in, D, D));
        
        l_tmp = deval(solution_l, t_in);
          
        matl_tmp = transpose(reshape(transpose(l_tmp), D, D));
        
        out_tmp=transpose(Q)*F_in+F_in*conj(Q)+transpose(R)*F_in*conj(R) ...
                + (1/(2*mass))*transpose(Rkin)*matl_tmp*conj(Rkin) ...
                + uu*(transpose(R)*transpose(R)*matl_tmp + matl_tmp*conj(R)*conj(R) )...
                + potential*transpose(R)*matl_tmp*conj(R) + interaction*transpose(R)*transpose(R)*matl_tmp*conj(R)*conj(R);
            %- energyDensity*matl_tmp;
        out = reshape(transpose(out_tmp), D^2, 1);

end

