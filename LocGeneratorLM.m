function [newLoc, count] = LocGeneratorLM(t, n, fieldDims, ...
            particleLocs, particleTypes, particleRs) 
%creates random loc for particle, ensuring it is not overlapping any
%current particles
    overlapping = 1;
    count = -1;
    while overlapping > 0
        newLoc = (2*rand(1,3)-1)*fieldDims(2,1);
        overlapping = CheckLocLM(newLoc, t, n, fieldDims, particleLocs, particleTypes, particleRs);
        count = count + 1;
    end
end
