clear;


fIn = sprintf('data/D=%d/cpxD=%dv=%dc=%d.mat',16, 16, 1, 100);

load(fIn);

p_min = 0; p_max = 4;

delta_p=0.5;

plot_x = p_min:delta_p:p_max;
plot_y = [];

T = kron(Q, eye(D,D))+ kron( eye(D,D), conj(Q)) + kron(R, conj(R));

proj = eye(D^2,D^2) - kron(r, l);


for (k = 1: numel(plot_x))
   fprintf('%d ', k);
   if (rem(k,10)==0)
        fprintf('\n');
   end
    
    Tp_minus = inv(  -T - kron(r, l)- i*plot_x(k)*eye(D^2, D^2) );
    Tp_plus = inv( -T- kron(r, l) + i*plot_x(k)*eye(D^2, D^2) );
    Tp_plus = proj*Tp_plus*transpose(proj);
    Tp_minus = proj*Tp_minus*transpose(proj);
    
    correlation_p = (1/occupationDensity)*( l*(kron(eye(D,D), conj(R))*Tp_plus*kron(R, eye(D,D)))*r ...
        + l*(kron(eye(D,D), conj(R))*Tp_minus*kron(R, eye(D,D)))*r );
    
    plot_y = [plot_y correlation_p];
    
end


figure; plot((1/occupationDensity)*plot_x, plot_y, '-xk');