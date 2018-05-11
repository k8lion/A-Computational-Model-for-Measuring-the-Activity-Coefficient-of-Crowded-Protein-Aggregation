function [rotStep] = RotStepGenerator(timeStep, t, n, particleRs, particleOrients) 
    k_B = 1.38064852e-23; %meters^2*kilograms/(seconds^2*kelvin)
    Temp = 298.15; %kelvin, body temperature
    eta = 8.91*10^-4; %Newton*seconds/meter^2 = kilograms/(meters*seconds)
    
    rotDiff = k_B*Temp/(8*eta*pi*(particleRs(n)*10^-10)^3); %1/seconds, Stokes-Einstein equation for spheres
    rotDiff = rotDiff*10^-12; %convert to 1/picoseconds

    rotStep = (2*rotDiff*timeStep)^(1/2)*randn(1,3);
end
