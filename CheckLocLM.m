function [overlapping] = CheckLocLM(tryLoc, t, n, fieldDims, ...
            particleLocs, particleTypes, particleRs)
    overlapping = 0;
    for m = 1:size(particleLocs, 1)
        if (m ~= n) & (DistanceBetween(particleLocs(m,:),tryLoc) < particleRs(n) + particleRs(m))
            overlapping = overlapping + 1;
        end
    end
end