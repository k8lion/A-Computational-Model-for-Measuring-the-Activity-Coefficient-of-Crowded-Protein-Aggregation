function [newLoc, stepDistSquared, count] = LocStepGenerator(timeStep, t, n, fieldDims, ...
            particleLocs, particleTypes, particleRs) 
%creates random step for particle, ensuring it is not overlapping any
%other particles
    k_B = 1.38064852e-23; %meters^2*kilograms/(seconds^2*kelvin)
    Temp = 298.15; %kelvin, body temperature
    eta = 8.91e-4; %Newton*seconds/meter^2 = kilograms/(meters*seconds)
    
    transDiff = k_B*Temp/(6*pi*eta*(particleRs(n)*10^-10)); %meters^2/seconds, Stokes-Einstein equation for spheres
    transDiff = transDiff*10^8; %convert to angstroms^2/picoseconds
    
    overlapping = 1;
    count = -1;
    while overlapping > 0   
        locStep = (2*transDiff*timeStep)^(1/2)*randn(1,3); %last term is time step in nanoseconds
        newLoc = particleLocs(n,:) + locStep;
        stepDistSquared = locStep(1)^2 + locStep(2)^2 + locStep(3)^2;
        for d = 1:3
           if newLoc(d) < fieldDims(1,d)
               newLoc(d) = newLoc(d) - 2*fieldDims(1,d);
           elseif newLoc(d) > fieldDims(2,d)
               newLoc(d) = newLoc(d) - 2*fieldDims(2,d);
           end
        end
        overlapping = CheckLocLM(newLoc, t, n, fieldDims, particleLocs, particleTypes, particleRs);
        count = count + 1;
    end
end






