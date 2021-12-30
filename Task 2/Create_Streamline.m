function [streamline, press_dist, x_crossings] = ...
Create_Streamline(initial_height_streamline, x_mesh, mach_initial,...
gamma, pressure_total, N, x_steps, intersections, gammas, Wave)

streamline = initial_height_streamline + 0 * x_mesh;
press_dist = (Pressure_From_Total_Pressure(mach_initial, gamma, pressure_total));
x_crossings = x_mesh(1);

if initial_height_streamline > 0
    end_reached = true;
else
    end_reached = false;
    [press_dist, x_crossings] = Center_Pressure_Distribution(press_dist, x_crossings);
end
prev_i = 1;

slope = 0;
while end_reached
    streamline_indices = Streamline_Crossing(streamline, prev_i, N, gammas);
    [prev_i, Gamma] = Find_Crossing(streamline_indices, prev_i, x_steps, N);
    if x_mesh(prev_i) < intersections.BCE_x(Gamma,Gamma)
            slope = tan(Wave.ABC_phi(Gamma));
            local_mach = Wave.ABC_mach(Gamma);
    elseif (x_mesh(prev_i) > intersections.BCE_x(Gamma,Gamma) && x_mesh(prev_i)...
            < intersections.BCE_x(Gamma,end))
        for i = 1:N
            if i == 1
                if (streamline(prev_i) < intersections.BCE_y(Gamma,i)...
                    && streamline(prev_i) > intersections.BCE_y(Gamma,end))
                    slope = tan((Wave.BCE_phi(Gamma,i) + Wave.BCE_phi(Gamma,end)) / 2);
                    local_mach = Wave.BCE_mach(Gamma,i);
                end
            else
                if (streamline(prev_i) < intersections.BCE_y(Gamma,i)...
                    && streamline(prev_i) > intersections.BCE_y(Gamma,i-1))
                    slope = tan((Wave.BCE_phi(Gamma,i) + Wave.BCE_phi(Gamma,i-1)) / 2);
                    local_mach = Wave.BCE_mach(Gamma,i);
                end
            end
        end
    elseif (x_mesh(prev_i) > intersections.BCE_x(Gamma,end) && x_mesh(prev_i)...
            < intersections.DFG_x(Gamma,Gamma))
            slope = tan(Wave.BCE_phi(Gamma,end));
            local_mach = Wave.BCE_mach(Gamma,end);
    elseif (x_mesh(prev_i) > intersections.DFG_x(Gamma,Gamma) && x_mesh(prev_i)...
            < intersections.DFG_x(Gamma,end))
        for i = 1:N
            if i == 1
                if (streamline(prev_i) < intersections.DFG_y(Gamma,i)...
                    && streamline(prev_i) > intersections.DFG_y(Gamma,end)...
                    || streamline(prev_i) > intersections.DFG_y(Gamma,i) &&...
                    streamline(prev_i) < intersections.DFG_y(Gamma,end))
                    slope = tan((Wave.DFG_phi(Gamma,i) + Wave.DFG_phi(Gamma,end))/2);
                    local_mach = Wave.DFG_mach(Gamma,i);
                end
            else
                if (streamline(prev_i) < intersections.DFG_y(Gamma,i)...
                   && streamline(prev_i) > intersections.DFG_y(Gamma,i-1)...
                   || streamline(prev_i) > intersections.DFG_y(Gamma,i)...
                   && streamline(prev_i) < intersections.DFG_y(Gamma,i-1))
                   slope = tan((Wave.DFG_phi(Gamma,i) + Wave.DFG_phi(Gamma,i-1))/2);
                   local_mach = Wave.DFG_mach(Gamma,i);
                end
            end
        end           
    elseif (x_mesh(prev_i) > intersections.DFG_x(Gamma,end)...
            && x_mesh(prev_i) < intersections.HIJ_x(Gamma,Gamma))
            slope = tan(Wave.DFG_phi(Gamma,end));
            local_mach = Wave.DFG_mach(Gamma,end);
    elseif (x_mesh(prev_i) > intersections.HIJ_x(Gamma,Gamma)...
            && x_mesh(prev_i) < intersections.HIJ_x(Gamma,end))
        for i = 1:N
            if i == 1
                if (streamline(prev_i) < intersections.HIJ_y(Gamma,i)...
                    && streamline(prev_i) > intersections.HIJ_y(Gamma,end))
                    slope = tan((Wave.HIJ_phi(Gamma,i) + Wave.HIJ_phi(Gamma,end))/2);
                    local_mach = Wave.HIJ_mach(Gamma,i);
                end
            else
                if (streamline(prev_i) < intersections.HIJ_y(Gamma,i)...
                    && streamline(prev_i) > intersections.HIJ_y(Gamma,i-1))
                    slope = tan((Wave.HIJ_phi(Gamma,i) + Wave.HIJ_phi(Gamma,i-1))/2);
                    local_mach = Wave.HIJ_mach(Gamma,i);
                end
            end
        end
    elseif x_mesh(prev_i) > intersections.HIJ_x(Gamma,end)
            slope = tan(Wave.HIJ_phi(Gamma,end));
            end_reached = false;
    end

    streamline(prev_i:end) = streamline(prev_i) + slope...
                             * (x_mesh(prev_i:end) - x_mesh(prev_i));
    
    x_crossings = [x_crossings, x_mesh(prev_i)];
    press_dist = [press_dist, Pressure_From_Total_Pressure(local_mach,...
                 gamma, pressure_total)];
end

end

