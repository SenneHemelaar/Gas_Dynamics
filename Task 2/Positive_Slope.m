function [positive_slope] = Positive_Slope(phi,mu)
%POSITIVE_SLOPE: Calculates positive slope from angles mu and phi in
%radians
positive_slope = tan(phi + mu);
end

