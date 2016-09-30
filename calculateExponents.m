%computes exponents
%run after structureFunction only for DSF and Gaussian

particle = false;

order = 0;

for k=1:numel(omega)
    if (omega(k) < excitation_energies(1))
        kOmegaHole = k;
    end
end
kOmegaHole = kOmegaHole+1;

[green_max, kOmegaParticle] = max(green);

%[curve, goodness] = fit(transpose(omega),transpose(green),'smoothingspline'); 
green_deconv = smooth(transpose(green)); green_deconv = transpose(green_deconv); dGreen = green_deconv;

if (particle == true)
    
    green_inv = zeros(1,numel(omega));
    
    Csquared_inv = 1./Csquared;
    
    for  m = 1: numel(omega)
   
        green_inv(m) = 0;
    
            prop_factor = 1;

        for k=1:D^2
            green_inv(m) = green_inv(m) ...
                +(1/fudge)*prop_factor*sqrt(pi/2)*Csquared_inv(k)*exp(-(omega(m) - excitation_energies(k))^2/(2*fudge^2));
        end
    end
    
    figure;
    plot(omega((kOmegaParticle-200):(kOmegaParticle+200))  , green_inv((kOmegaParticle-200):(kOmegaParticle+200)), 'x' );

elseif (particle == false)
    
    %HOLE log vs. log
    %figure;
   % hold on;
    %plot(omega/((pi*real(occupationDensity))^2), log(green_deconv) , 'x');
    % plot(omega(1:(end-1000))/((pi*occupationDensity)^2), green_deconv(1:(end-1000)) );
     log_epsilon = exp(-12);
     lowest_excitation = real(excitation_energies(1));
   % plot(log(omega((kOmegaHole+20):(kOmegaHole+50))-lowest_excitation), log(abs(green((kOmegaHole+20):(kOmegaHole+50))-green(kOmegaHole))), 'x' );
   % hold off;
    
    %PARTICLE log vs. log
    figure;
   plot(log(( omega((kOmegaParticle):(kOmegaParticle+50))-omega(kOmegaParticle) ) ), log( green((kOmegaParticle-50):(kOmegaParticle))  - green(kOmegaParticle)), 'x' );
    
end


    %deconvolution?
    %for k=1:order
    %    dy = diff(dGreen);
    %    dGreen = dy./delta_omega; dGreen = [dGreen dGreen(end) ] ;
    %    dGreen(1:order) = 0; dGreen((end-order):end) = 0;
    %    dGreen = smooth(transpose(dGreen), 'rloess'); dGreen = transpose(dGreen);
    
    %    dy = diff(dGreen);
    %    dGreen = dy./delta_omega; dGreen = [dGreen dGreen(end) ] ;
    %    dGreen(1:2*order) = 0; dGreen((end-2*order):end) = 0;
    %    dGreen = smooth(transpose(dGreen), 'rloess'); dGreen = transpose(dGreen);
            
    %    green_deconv = green_deconv + ( ((-1)^k)*( fudge^(2*k) )/( (2^k)*factorial(k) ) )*dGreen;

    %    green_deconv(1:2*order) = 0; green_deconv((end-2*order):end) = 0;
    %    green_deconv = smooth(transpose(green_deconv), 'rloess'); green_deconv = transpose(green_deconv);
    %end