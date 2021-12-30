function [press_dist, x_crossings] = Center_Pressure_Distribution(press_dist,x_crossings, Wave_BCE_mach, Wave_HIJ_mach, gamma, p_total)
%CENTER_PRESSURE_DISTRIBUTION Summary of this function goes here
for i = 1:N
    x_pos       = x_intersections_BCE(i,i);
    press_dist  = [ [press_dist],  [Pressure_From_Total_Pressure(Wave_BCE_mach(i,i), gamma, p_total)] ];
    x_crossings = [ [x_crossings] , [x_pos] ];
end

for i = 1:N
    x_pos       = x_intersections_HIJ(i,i);
    press_dist  = [ [press_dist] , [Pressure_From_Total_Pressure(Wave_HIJ_mach(i,i), gamma, p_total)] ];
    x_crossings = [ [x_crossings] , [x_pos] ];
end

end

