function [newLoc] = StepGenerator(t, n, fieldDims, ...
            particleLocs, particleTypes, particleRs) 
%creates random step for particle, ensuring it is not overlapping any
%other particles
    k_B = 1.38064852e-23; %meters^2*kilograms/(seconds^2*kelvin)
    Temp = 298.15; %kelvin, body temperature
    eta = 8.91e-4; %Newton*seconds/meter^2 = kilograms/(meters*seconds)
    
    transDiff = k_B*Temp/(6*pi*eta*(particleRs(n)*10e-10)); %meters^2/seconds, Stokes-Einstein equation for spheres
    transDiff = transDiff*10e11; %convert to angstroms^2/nanoseconds
    
    overlapping = 1;
    while overlapping > 0   
        locStep = 2*transDiff*randn(1,3) 
        newLoc = particleLocs(n,:,t-1) + locStep;
        for d = 1:3
           if newLoc(d) < fieldDims(1,d)
               newLoc(d) = newLoc(d) - 2*fieldDims(1,d);
           elseif newLoc(d) > fieldDims(2,d)
               newLoc(d) = newLoc(d) - 2*fieldDims(2,d);
           end
        end
        overlapping = CheckLoc(newLoc, t, n, fieldDims, particleLocs, particleTypes, particleRs);
    end
end






