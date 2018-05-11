function [dist] = DistanceBetween(loc1, loc2)
%returns distance between 2 locations
    dist = sqrt((loc1(1)-loc2(1))^2+(loc1(2)-loc2(2))^2+(loc1(3)-loc2(3))^2);
end
