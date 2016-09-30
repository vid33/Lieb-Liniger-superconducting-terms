 function out=T_action(in, R, Q, D)
        in = transpose(reshape(in, D, D));
        out=Q*in+in*Q'+R*in*R';
        out = reshape(transpose(out), D^2, 1);
 end

