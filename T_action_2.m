function out=T_action_2(t_in, x_in, R, Q, D)
        x_in = transpose(reshape(x_in, D, D));
        out_tmp=transpose(Q)*x_in+x_in*conj(Q)+transpose(R)*x_in*conj(R);
        out = reshape(transpose(out_tmp), D^2, 1);
        
        %out = transpose(out);
end

