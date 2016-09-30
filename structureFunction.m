clear;

D=16;
potential =  10;
interaction = 0;
uu = -5;

DSF = false;
Gaussian = true; %otherwise convolute with Lorentzian
testNorm  = false; %compare sum of squared form factors with norm; slow for high D

%omega_start = 1.65; omega_end = 2;
omega_start=0; omega_end=505;
delta_omega=0.005;
omega = omega_start:delta_omega:omega_end;

fudge = 0.1; %width of lorentzian/gaussian
green = zeros(1,numel(omega));

fIn_density = sprintf('data_LLuu/D=%d/D=%duu=%dv=%dc=%d.mat', D, D, uu, potential, interaction);

load(fIn_density, 'occupationDensity');

p = 0*occupationDensity;

if (DSF == true)
    fIn = sprintf('data_LLuu/D=%d/DSF_form_factors/DSF_form_factorsD=%duu=%dv=%dc=%dp=%d.mat', D, D, uu, potential, interaction, p);
elseif (DSF == false)
    fIn = sprintf('data_LLuu/D=%d/form_factors/form_factorsD=%duu=%dv=%dc=%dp=%d.mat', D, D, uu, potential, interaction, p);
end  
    
load(fIn);

Csquared = zeros(1, D^2);
for k=1:D^2
    Csquared(k) = C{k}*conj(C{k});
end

if (DSF == false)
    C1Psquared = zeros(1, D^2);
    for k=1:D^2
        C1Psquared(k) = C1P{k}*conj(C1P{k});
    end
end
 
for  m = 1: numel(omega)
   
    green(m) = 0;
    
    if (DSF == true)
        prop_factor = 1;
    elseif (DSF == false)
        prop_factor = (1/pi);
    end
    % hold all;
   % plot((excitation_energies./(2*pi*occupationDensity^2)), Csquared, 'x');
   % hold off;    
    %figure('Name', figureName); plot(omega, green);
    
    for k=1:D^2
        
        if (Gaussian == true)
            if (omega(m) > 0)
                green(m) = green(m) ...
                    +(1/fudge)*prop_factor*sqrt(pi/2)*Csquared(k)*exp(-(omega(m) - excitation_energies(k))^2/(2*fudge^2));
            elseif (omega(m) < 0)
                green(m) = green(m) ...
                    +(1/fudge)*prop_factor*sqrt(pi/2)*C1Psquared(k)*exp(-(-omega(m) - excitation_energies(k))^2/(2*fudge^2));
            end
        else %Lorentzian
            if (omega(m) > 0)
                green(m) = green(m) ...
                    + prop_factor*Csquared(k)*fudge/( (omega(m)-excitation_energies(k))^2 + fudge^2);
            elseif (omega(m) < 0)
                green(m) = green(m) ...
                    + prop_factor*C1Psquared(k)*fudge/( (-omega(m)-excitation_energies(k))^2 + fudge^2);
            end
        end
    end
end

if (DSF == true)
    fsum_rule_lhs = 0;
    for k = 1: D^2
        fsum_rule_lhs = fsum_rule_lhs + excitation_energies(k)*C{k}*conj(C{k});
    end
end

figureName = sprintf('D=%d, c_eff=%g, p_eff=%g, fudge=%g, domega=%g', D, interaction_effective, p/occupationDensity, fudge, delta_omega);

if (DSF == true)
    figure('Name', figureName); 
    %hold all
    plot(omega, green);    
    %plot(omega/((pi*occupationDensity)^2), green);    
    %hold off;
    fprintf('f-sum lhs is %d, rhs is %d \n', fsum_rule_lhs, occupationDensity*p^2);
elseif (DSF == false)
    figure('Name', figureName); 
    %hold all;
    log_epsilon = exp(-12);
    plot(omega/((pi*occupationDensity)^2), log(real(green)+log_epsilon) );
    %plot(omega/((pi*occupationDensity)^2), green );
    %hold off;
end


%Plot form factors
if (DSF == true)   
    kmax = 0; kOmegaHole=0;
    for k=1:D^2
        if (excitation_energies(k)  <= omega_end)
            kmax = k;
        end
    end
    
    
    %hold all;
    %stem(excitation_energies(1:kmax)/((pi*occupationDensity)^2), Csquared(1:kmax), 'x' );
    %stem(excitation_energies(1:kmax), Csquared(1:kmax), 'x' );
    %hold off;
end
 
    
%%%%%%%%%%%%%%%%%%%% STUFF THAT CHECKS NORM vs. SUM OF SQUARED FORM
%%%%%%%%%%%%%%%%%%%% FACTORS; EXPENSIVE

if (testNorm == true)
    
    proj = eye(D^2,D^2) - kron(r, l);

    Tplus = - kron(Q, eye(D,D)) - kron( eye(D,D), conj(Q)) - kron(R, conj(R)) - kron(r,l) + 1i*p*eye(D^2, D^2);
    Tminus = - kron(Q, eye(D,D)) - kron( eye(D,D), conj(Q)) - kron(R, conj(R)) - kron(r,l) - 1i*p*eye(D^2, D^2);
    if (DSF == true)
        DSF_norm  = (l*(kron(R, conj(R))*proj)/Tminus)*proj*kron(R, conj(R))*r ...
            +(l*(kron(R, conj(R))*proj)/Tplus)*proj*kron(R, conj(R))*r + l*kron(R, conj(R))*r;
        fprintf('Norm of psi dagger(psi) %6.14f vs. sum form_factors^2  %6.14f \n', DSF_norm, sum(Csquared));
    elseif (DSF == false)
        spectralParticle_norm  = (l*(kron(R, eye(D,D))*proj)/Tminus)*proj*kron(eye(D,D), conj(R))*r ...
            +(l*(kron(eye(D,D), conj(R))*proj)/Tplus)*proj*kron(R, eye(D,D))*r + l*r;
        spectralHole_norm  = (l*(kron(R, eye(D,D))*proj)/Tplus)*proj*kron(eye(D,D), conj(R))*r ...
            +(l*(kron(eye(D,D), conj(R))*proj)/Tminus)*proj*kron(R, eye(D,D))*r;
         fprintf('Norm of dagger(psi) %6.14f vs. sum form_factors^2 %6.14f \n', spectralParticle_norm, sum(Csquared));
         fprintf('Norm of psi %6.14f vs. sum form_factors^2 %6.14f \n', spectralHole_norm, sum(C1Psquared));
    end
   
end