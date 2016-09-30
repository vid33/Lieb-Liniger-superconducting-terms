function out=T_action_left_const(t_in, x_in, R, Q, D)
        x_in = transpose(reshape(x_in, D, D));
        out_tmp=-(Q*x_in + x_in*Q' + R*x_in*R');
        
        out_tmp = 0.5*(out_tmp+ out_tmp');
        
        out = reshape(transpose(out_tmp), D^2, 1);
        
        %out = transpose(out);
end

