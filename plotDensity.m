function [ output ] = plotDensity( density, N)

    for k = 1:(N)
        plotData(k) = real(density{k});
    end
    output = 0;
    
    plot(plotData);
    
end

