function [newPoint] = Rotate(oldPoint, rotStep) 
    xRot = [1 0 0;0 cos(rotStep(1)) -sin(rotStep(1));0 sin(rotStep(1)) cos(rotStep(1))]; %rotation around x axis
    yRot = [cos(rotStep(2)) 0 sin(rotStep(2));0 1 0;-sin(rotStep(2)) 0 cos(rotStep(2))]; %rotation around y axis
    zRot = [cos(rotStep(3)) -sin(rotStep(3)) 0;sin(rotStep(3)) cos(rotStep(3)) 0;0 0 1]; %rotation around z axis
    
    newPoint = oldPoint*((zRot*yRot*xRot)');
end
