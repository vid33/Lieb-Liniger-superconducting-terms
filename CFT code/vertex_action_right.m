function out=vertex_action(t_in, x_in, beta, Lambda, R, D)
        x_in = transpose(reshape(x_in, D, D));
        out_tmp=(-beta/sqrt(2*Lambda)  )*1i*( R*x_in+x_in*R' );
        %+transpose(R)*x_in*conj(R);
        out = reshape(transpose(out_tmp), D^2, 1);
        
        %out = transpose(out);
end

