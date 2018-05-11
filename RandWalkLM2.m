%Kaitlin Maile
%Wood Group
%RandWalk

function [timeSteps, particleLocs, particleTypes, particleRs, ...
    particleOrients, volFracA, totalTime, transDiff, effectiveDiff, ...
    bondPointsR, numBonds, timesToBond, conc, coll, collCount] = RandWalkLM2(reactantR, crowderR, ...
    fieldRatio, numReactants, numCrowders, totalSteps, maxTotalSteps, initBond, det)

    tstart = tic;

    k_B = 1.38064852e-23; %meters^2*kilograms/(seconds^2*kelvin)
    temperature = 298.15; %kelvin, body temperature
    eta = 8.91e-4; %Newton*seconds/meter^2 = kilograms/(meters*seconds)

    fieldDims = fieldRatio*[-reactantR, -reactantR, -reactantR;...
                             reactantR,  reactantR,  reactantR]; %initiates 
                         %field: cubic, based on center particle's radius 
    conc = -numReactants/(fieldDims(1,1)*2)^3; %calculated concentration                   
    totalTime = 0; %elapsed time in picoseconds
    timesToBond = [0,0,0,0]; %time to reach N=1,2,3,4 states
    sumSquareDisp = zeros(numReactants,1); %tracks sum squared displacement for each reactant
    sumDisp = zeros(numReactants,1); %tracks sum displacement for each reactant
    allEffArcAngles = zeros(numReactants,1); %tracks effective arcs
    allActArcAngles = zeros(numReactants,3,1); %tracks actual effective arcs
    sumSquareRot = zeros(numReactants,3); %tracks sum squared rot for each reactant
    timeStep = .2*det; %ps, first/default timestep
    timeSteps = timeStep; %1 x 1 x totalSteps, pages = time step, in ps                   
    particleLocs = zeros(1,3); %rows = particles, columns = dimensions (range def. by fieldDims)
    particleOrients = zeros(numReactants,3); %rows = particles, columns = dimensions (range -pi to pi)
    particleTypes = zeros(numReactants+numCrowders,1); %rows = particles; 1 = reactive (reactant), 0 = inert (crowder)
    particleRs = zeros(numReactants+numCrowders,1); %rows = particles, value is that particle's radius
    bondR = reactantR*.45; %for proper placement of bond points on reactant surface
    bondPointsR = zeros(4,3,numReactants); %bond point, dimension, reactant number
    bondDists = zeros(4,1,numReactants);
    lockPot = -4.2; %kcal/mol; locking potential energy
    dissocProb = exp(lockPot/(.001986*temperature)); %probability of dissociation
    numBonds = zeros(numReactants, totalSteps); %tracks number of bonds to center reactant for each time step
    detailCoeff = det; %higher -> shorter timestep, more detailed but more intensive
    transDiff = k_B*temperature/(6*pi*eta*(reactantR*10^-10)); %meters^2/seconds, Stokes-Einstein equation for spheres
    transDiff = transDiff*10^8; %convert to angstroms^2/picoseconds
    actRot = 10^-12*k_B*temperature/(8*eta*pi*(reactantR*10^-10)^3); %1/picoseconds, Stokes-Einstein equation for spheres
    coll = [0];
    for n = 1:numReactants+numCrowders %for each particle
        %first, fill in relevant initial data: type, radius, orientation
        if n > numReactants %if a crowder
            particleTypes(n) = 0;
            particleRs(n) = crowderR;
        elseif n == 2 %if second reactant
            particleTypes(n) = 1;
            particleRs(n) = reactantR;
            
            if initBond == 4
                particleOrients(n,:) = [-pi,0,pi/2]; %N=4 fully bonded
            elseif initBond == 2
                particleOrients(n,:) = [-pi,pi/2,pi/2]; %N=2
            elseif initBond == 1
                particleOrients(n,:) = [-3*pi/4,pi/4,pi/2]; %N=1
            elseif initBond == 0
                particleOrients(n,:) = (2*rand(1,3)-1)*pi; %near field regime
            else
                particleOrients(n,:) = (2*rand(1,3)-1)*pi; %random orient
            end
        else %if any other reactant
            particleTypes(n) = 1;
            particleRs(n) = reactantR;
            particleOrients(n,:) = (2*rand(1,3)-1)*pi; %random orient
        end
        if initBond > 0
            particleOrients(1,:) = [0,0,0]; %fix initial orient of center 
            %reactant if needed for second reactant placement purposes
        end
        
        %next, determine initial location for this particle
        if n == 1
            nthParticleLoc = [0,0,0];
        elseif n == 2
            if initBond == 4
                nthParticleLoc = [0,0,2.05*reactantR]; %N=4 fully bonded
            elseif initBond == 2
                nthParticleLoc = [0,1.43*reactantR,1.43*reactantR]; %N=2
            elseif initBond == 1
                nthParticleLoc = [.85*reactantR,.85*reactantR,2.05*reactantR]; %N=1
            elseif initBond == 0
                nthParticleLoc = [0,0,42]; %near field regime  
            else
                [nthParticleLoc, ~] = LocGeneratorLM(1, n, fieldDims, ...
                    particleLocs, particleTypes, particleRs); %random loc
            end
        else 
            [nthParticleLoc, ~] = LocGeneratorLM(1, n, fieldDims, ...
                particleLocs, particleTypes, particleRs);
        end
        particleLocs(n,:) = nthParticleLoc;
    end
    
    bondPointsR(1,:,1) = [bondR,bondR,reactantR]; %relative coords of bond points for orient of 0
    bondPointsR(2,:,1) = [bondR,-bondR,reactantR];
    bondPointsR(3,:,1) = [-bondR,-bondR,reactantR];
    bondPointsR(4,:,1) = [-bondR,bondR,reactantR];

    numBond = zeros(numReactants,1); %temporary to keep track for current time step then update numBonds
    
    for n = 2:numReactants
        bondPointsR(1,:,n) = [bondR,bondR,reactantR]; %relative coords of bond points for orient of 0
        bondPointsR(2,:,n) = [-bondR,bondR,reactantR];
        bondPointsR(3,:,n) = [-bondR,-bondR,reactantR];
        bondPointsR(4,:,n) = [bondR,-bondR,reactantR];
        for b = 1:4 %adjust bond points to initial orientation
            bondPointsR(b,:,n) = Rotate(bondPointsR(b,:,n), particleOrients(n,:));
            bondDists(b,1,n) = DistanceBetween(particleLocs(n,:)+bondPointsR(b,:,n),particleLocs(1,:)+bondPointsR(b,:,1));
                %find distance between center reactant BPs and this
                %reactant's BPs
            if bondDists(b,1,n) < 2  
                numBond(n) = numBond(n) + 1;
            end
        end
    end
    
    numBonds(:,1) = numBond(:);
    
    %calc vol frac by setting up random grid of points and determining how
    %many are in a particle
    A = (2*rand(10^4,3)-1)*fieldDims(1,1);
    volFracA = sum(sum(pdist2(A,particleLocs(1:numReactants,:))<=reactantR))/(size(A,1))...
               +sum(sum(pdist2(A,particleLocs(numReactants+1:numReactants+numCrowders,:))<=crowderR))/(size(A,1));
    %volFracA = numReactants*pi*(2*reactantR)^3/(6*(2*fieldDims(2,1)+2*reactantR)^3); %Yeth method
    
    collCount = 0;
    collBool = 0;
    t = 1;
    %run rest of steps
    while ((sum(numBonds(:,t) < 4) || t < totalSteps) && (t < maxTotalSteps) && (toc(tstart) < 90*60*60))
        %use miminum distance to determine variable time step
        minDist = 2*fieldDims(2,1);
        for p2 = 2:numReactants+numCrowders
            currDist = DistanceBetween(particleLocs(1,:),particleLocs(p2,:));
            if currDist < minDist
                minDist = currDist;
            end
        end
        if minDist < 2 + 2*reactantR && all(numBonds(:) <  1)
            collBool = 1;
            if coll(t) == 0
                collCount = collCount + 1;
            end
        elseif minDist > 4 + 2*reactantR && all(numBonds(:) <  1)
            collBool = 0;
        end
            
        coll(t+1) = collBool;
                
        t = t + 1;
        if all(numBonds(:) <  1)
            %determine current minimum distance between particles
            minDist = 2*fieldDims(2,1);
            for p1 = 1:numReactants+numCrowders-1
                for p2 = p1+1:numReactants+numCrowders
                    currDist = DistanceBetween(particleLocs(p1,:),particleLocs(p2,:));
                    if currDist < minDist
                        minDist = currDist;
                    end
                end
            end
            timeSteps(1,1,t) = (minDist/detailCoeff)^2/(12*transDiff);
            %use miminum distance to determine variable time step
        else
            timeSteps(1,1,t) = .2;
        end
        totalTime = totalTime + timeSteps(1,1,t);
        currTime(t) = totalTime;
        
        %first reactant does not move
        newOrientStep = RotStepGenerator(timeSteps(1,1,t), t, 1, particleRs, particleOrients);
        particleOrients(1,:) = newOrientStep + particleOrients(1,:);
        for b = 1:4
            bondPointsR(b,:,1) = Rotate(bondPointsR(b,:,1), newOrientStep);
        end
        allEffArcAngles(1,t) = reactantR*acos(dot([1,0,0], Rotate([1,0,0], newOrientStep)));
        allActArcAngles(1,:,t) = reactantR*newOrientStep;
        sumSquareRot(1,:) = sumSquareRot(1,:) + newOrientStep.^2;
        for n = 2:numReactants+numCrowders %for rest of particles
            %find and assign next location
            if particleTypes(n) == 1 %if particle is reactant
                stepAccepted = false;
                while stepAccepted == false
                    numBondn = 0;
                    [newLoc, dispSqr, ~] = LocStepGeneratorLM(timeSteps(1,1,t), t, n, fieldDims, ...
                        particleLocs, particleTypes, particleRs);
                    %rejSteps(t) = counter;
                
                    %find and assign next orientation
                    newOrientStep = RotStepGenerator(timeSteps(1,1,t), t, n, particleRs, particleOrients);
                    tempBondPointsR = zeros(4,3);
                    tempBondDists = zeros(4,1);
                    
                    for b = 1:4
                        tempBondPointsR(b,:) = Rotate(bondPointsR(b,:,n), newOrientStep);
                        tempBondDists(b) = DistanceBetween(newLoc+tempBondPointsR(b,:),particleLocs(1,:)+bondPointsR(b,:,1));
                        if tempBondDists(b) < 2 
                            numBondn = numBondn + 1;
                        end
                    end
                
                    numBonds(n,t) = numBondn;
                
                    if numBonds(n,t) < 2 || (rand < dissocProb || numBonds(n,t) - numBonds(n,t-1) >= 0)     
                        for b = 1:4
                            bondPointsR(b,:,n) = tempBondPointsR(b,:);
                            bondDists(b,1,n) = tempBondDists(b);
                        end
                        particleOrients(n,:) = particleOrients(n,:) + newOrientStep;
                        particleLocs(n,:) = newLoc;
                    
                        %add on step to total distance traveled
                        allEffArcAngles(n,t) = reactantR*acos(dot([1,0,0], Rotate([1,0,0], newOrientStep)));
                        allActArcAngles(n,:,t) = reactantR*newOrientStep;
                        sumSquareDisp(n) = sumSquareDisp(n) + dispSqr;
                        sumDisp(n) = sumDisp(n) + dispSqr^(1/2);
                        sumSquareRot(n,:) = sumSquareRot(n,:) + newOrientStep.^2;
                        stepAccepted = true;
                    end
                end
                
            else
                [newLoc, ~, ~] = LocStepGeneratorLM(timeSteps(1,1,t), t, n, fieldDims, ...
                    particleLocs, particleTypes, particleRs);
                particleLocs(n,:) = newLoc;
            end
        end
    end
    effectiveDiff = sumSquareDisp(:)/(6*sum(timeSteps));
    effectiveDiff(1) = transDiff;
    effectiveRot = sumSquareRot(:,:)/(2*sum(timeSteps));
    for b = 1:4
        [~,j] = find(sum(numBonds) == b,1);
        if j > 1
            timesToBond(b) = currTime(j);
        else
            timesToBond(b) = totalTime;
        end
    end
    timesToBond(5) = totalTime;
    %sumAngles = sumAngles./(totalTime*.001)
end
 