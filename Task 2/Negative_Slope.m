function [negative_slope] = Negative_Slope(phi,mu)
%POSITIVE_SLOPE: Calculates negative slope from angles mu and phi in
%radians
negative_slope = tan(phi - mu);
end

