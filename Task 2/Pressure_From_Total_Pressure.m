function [p] = Pressure_From_Total_Pressure(mach, gamma, p_total)
%PRESSURE_FROM_TOTAL_PRESSURE
p = (1 + (gamma-1)/2 * mach^2)^(-gamma/(gamma-1)) * p_total;
end

