%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       Gas Dynamics: Task 2                          %%%
%%%                      Main processing script                         %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%                           October 2021                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%%%                            CONSTANTS                                %%%
%%%=====================================================================%%%
gamma = 7/5;   % Must also change in function 'Prandtl_Meyer_From_Mach.m'

%%%                          INITIAL VALUES                             %%%
%%%=====================================================================%%%
mach_initial = 1.95;
phi_initial = 0;
pressure_ambient = 1225;
pressure_initial = 2*pressure_ambient;
pressure_total = 1;

%%%                        VARIABLE PARAMETERS                          %%%
%%%=====================================================================%%%
N = 15;
h = 1;
x_range = 14;
x_steps = 10000;
x_mesh  = linspace(0, x_range, x_steps);
error_tolerance = 0.00001;

%%%                      PRILIMINARY CALCULATIONS                       %%%
%%%=====================================================================%%%

mach_B        = sqrt(2 / (gamma-1) * (2^((gamma - 1) / gamma) *...
                (1 + (gamma-1) / 2 * mach_initial^2) -1));
phi_B         = Forward_Phi_Gamma_Plus(Prandtl_Meyer_From_Mach(mach_initial),...
                Prandtl_Meyer_From_Mach(mach_B), phi_initial);
alpha_initial = Alpha(Mu_From_Mach(mach_initial),phi_initial);
alpha_B       = Alpha(Mu_From_Mach(mach_B),phi_B);
alpha_list    = linspace(alpha_initial,alpha_B,N);
alpha_steps   = linspace(1,length(alpha_list),length(alpha_list));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Setting up zero arrays and matrices for all regions         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave.ABC_mach  = zeros(N,1);
Wave.ABC_nu    = zeros(N,1);
Wave.ABC_phi   = zeros(N,1);

Wave.BCE_mach         = zeros(N);
Wave.BCE_nu           = zeros(N);
Wave.BCE_phi          = zeros(N);
Wave.BCE_alpha        = zeros(N);
Wave.BCE_mu           = zeros(N);
intersections.BCE_x   = zeros(N);
intersections.BCE_y   = zeros(N);

Wave.DFG_mach         = zeros(N);
Wave.DFG_nu           = zeros(N);
Wave.DFG_phi          = zeros(N);
Wave.DFG_alpha        = zeros(N);
Wave.DFG_mu           = zeros(N);
intersections.DFG_x   = zeros(N);
intersections.DFG_y   = zeros(N);

Wave.HIJ_mach         = zeros(N);
Wave.HIJ_nu           = zeros(N);
Wave.HIJ_phi          = zeros(N);
Wave.HIJ_alpha        = zeros(N);
Wave.HIJ_mu           = zeros(N);
intersections.HIJ_x   = zeros(N);
intersections.HIJ_y   = zeros(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Finding relevant variables* for all regions              %%%
%%%                  *: Mach, nu, phi, alpha, mu                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave.ABC_mach(1)  = mach_initial;
Wave.ABC_nu(1)    = Prandtl_Meyer_From_Mach(mach_initial);
Wave.ABC_phi(1)   = phi_initial;

%%%                         REGION ABC, BCE                             %%%
%%%=====================================================================%%%
for j = 2:N
    [Wave.ABC_mach(j), Wave.ABC_phi(j), Wave.ABC_nu(j)] =...
    Bisection_Method(...
    Wave.ABC_mach(j-1),...
    mach_B,...
    error_tolerance,...
    alpha_list(j),...
    Wave.ABC_phi(j-1),...
    Prandtl_Meyer_From_Mach(Wave.ABC_mach(j-1)));
    Wave.BCE_nu(1,1) = Wave.ABC_nu(1) + Wave.ABC_phi(1);
    for i = 1:j
        if i == j
            Wave.BCE_nu(i,j) = Wave.ABC_nu(j) + Wave.ABC_phi(j);
        elseif i == 1
            Wave.BCE_nu(i,j) =  0.5*(Wave.BCE_nu(i,j-1) + Wave.ABC_nu(j))+...
                                0.5*(Wave.ABC_phi(j) - Wave.BCE_phi(i,j-1));
            Wave.BCE_phi(i,j) = 0.5*(Wave.BCE_phi(i,j-1) + Wave.ABC_phi(j))+...
                                0.5*(Wave.ABC_nu(j) - Wave.BCE_nu(i,j-1));
        else
            Wave.BCE_nu(i,j) =  0.5*(Wave.BCE_nu(i,j-1) + Wave.BCE_nu(i-1,j))+...
                                0.5*(Wave.BCE_phi(i-1,j) - Wave.BCE_phi(i,j-1));
            Wave.BCE_phi(i,j) = 0.5*(Wave.BCE_phi(i,j-1) + Wave.BCE_phi(i-1,j))+...
                                0.5*(Wave.BCE_nu(i-1,j) - Wave.BCE_nu(i,j-1));
        end
        [Wave.BCE_mu(1,1), Wave.BCE_mach(1,1)] = Mu_From_Nu_Bisection(...
                                                 Wave.BCE_nu(1,1), mach_initial-0.007,... 
                                                 3, error_tolerance);
        [Wave.BCE_mu(i,j), Wave.BCE_mach(i,j)] = Mu_From_Nu_Bisection(...
                                                 Wave.BCE_nu(i,j), mach_initial,...
                                                 20, error_tolerance);
        Wave.BCE_alpha(1,1) = Positive_Slope(Wave.BCE_mu(1,1), Wave.BCE_phi(1,1));
        Wave.BCE_alpha(i,j) = Positive_Slope(Wave.BCE_mu(i,j), Wave.BCE_phi(i,j));
    end
end

%%%                           REGION DFG                                %%%
%%%=====================================================================%%%
for j = 2:N 
    for i = 1:j
        Wave.DFG_nu(1,1)  = Prandtl_Meyer_From_Mach(mach_B);
        Wave.DFG_phi(1,1) = Wave.BCE_phi(1,end) - Wave.BCE_nu(1,end) + Wave.DFG_nu(1,1);
        if i == j
            Wave.DFG_nu(i,j)  = Prandtl_Meyer_From_Mach(mach_B);
            Wave.DFG_phi(i,j) = Wave.DFG_phi(i-1,j) - Wave.DFG_nu(i-1,j) + Wave.DFG_nu(i,j);
        elseif i == 1
            Wave.DFG_nu(i,j)  = 0.5*(Wave.BCE_nu(j,end)  + Wave.DFG_nu(i,j-1))+...
                                0.5*(Wave.DFG_phi(i,j-1) - Wave.BCE_phi(j,end));
            Wave.DFG_phi(i,j) = 0.5*(Wave.BCE_phi(j,end) + Wave.DFG_phi(i,j-1))+...
                                0.5*(Wave.DFG_nu(i,j-1)  - Wave.BCE_nu(j,end));
        else
            Wave.DFG_nu(i,j)  = 0.5*(Wave.DFG_nu(i,j-1)  + Wave.DFG_nu(i-1,j))+...
                                0.5*(Wave.DFG_phi(i,j-1) - Wave.DFG_phi(i-1,j));
            Wave.DFG_phi(i,j) = 0.5*(Wave.DFG_phi(i,j-1) + Wave.DFG_phi(i-1,j))+...
                                0.5*(Wave.DFG_nu(i,j-1)  - Wave.DFG_nu(i-1,j));
        end
        [Wave.DFG_mu(1,1), Wave.DFG_mach(1,1)] = Mu_From_Nu_Bisection(...
                                                 Wave.DFG_nu(1,1), mach_initial,...
                                                 20, error_tolerance);
        [Wave.DFG_mu(i,j), Wave.DFG_mach(i,j)] = Mu_From_Nu_Bisection(...
                                                 Wave.DFG_nu(i,j), mach_initial,...
                                                 20, error_tolerance);
        Wave.DFG_alpha(1,1) = -Negative_Slope(Wave.DFG_mu(1,1), Wave.DFG_phi(1,1));
        Wave.DFG_alpha(i,j) = -Negative_Slope(Wave.DFG_mu(i,j), Wave.DFG_phi(i,j));
    end
end

%%%                           REGION HIJ                                %%%
%%%=====================================================================%%%
for j = 2:N 
    for i = 1:j
        Wave.HIJ_nu(1,1) = Wave.DFG_nu(1,end) + Wave.DFG_phi(1,end);
        if i == j
            Wave.HIJ_nu(i,j)  = Wave.DFG_nu(j,end) + Wave.DFG_phi(j,end);
        elseif i == 1
            Wave.HIJ_nu(i,j)  = 0.5*(Wave.HIJ_nu(i,j-1)  + Wave.DFG_nu(j,end))+...
                                0.5*(Wave.DFG_phi(j,end) - Wave.HIJ_phi(i,j-1));
            Wave.HIJ_phi(i,j) = 0.5*(Wave.HIJ_phi(i,j-1) + Wave.DFG_phi(j,end))+...
                                0.5*(Wave.DFG_nu(j,end)  - Wave.HIJ_nu(i,j-1));
        else
            Wave.HIJ_nu(i,j)  = 0.5*(Wave.HIJ_nu(i,j-1)  + Wave.HIJ_nu(i-1,j))+...
                                0.5*(Wave.HIJ_phi(i-1,j) - Wave.HIJ_phi(i,j-1));
            Wave.HIJ_phi(i,j) = 0.5*(Wave.HIJ_phi(i,j-1) + Wave.HIJ_phi(i-1,j))+...
                                0.5*(Wave.HIJ_nu(i-1,j)  - Wave.HIJ_nu(i,j-1));
        end
        [Wave.HIJ_mu(1,1), Wave.HIJ_mach(1,1)] = Mu_From_Nu_Bisection(...
                                                 Wave.HIJ_nu(1,1), mach_initial-0.007,...
                                                 20, error_tolerance);
        [Wave.HIJ_mu(i,j), Wave.HIJ_mach(i,j)] = Mu_From_Nu_Bisection(...
                                                 Wave.HIJ_nu(i,j), mach_initial-0.007,...
                                                 20, error_tolerance);
        Wave.HIJ_alpha(1,1) = Positive_Slope(Wave.HIJ_mu(1,1), Wave.HIJ_phi(1,1));
        Wave.HIJ_alpha(i,j) = Positive_Slope(Wave.HIJ_mu(i,j), Wave.HIJ_phi(i,j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Finding intersection coordinates and create plotting        %%%
%%%                   arrays for regions BCE, DFG, HIJ                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%                        REGION ABC, BCE                              %%%
%%%=====================================================================%%%
gammas = zeros(N,length(x_mesh));
for i = 1:N
    for j = i:N
        positive_slope = 0;
        negative_slope = 0;
        if i == 1
            negative_slope = -Negative_Slope(Mu_From_Mach(Wave.ABC_mach(j)),...
                             Wave.ABC_phi(j));
            incoming_characteristic = h + negative_slope * x_mesh;
            if j == 1
                reflected_characteristic = 0 * x_mesh;
            else
                positive_slope = Positive_Slope(Wave.BCE_mu(i,j), Wave.BCE_phi(i,j));
                diff_positive_slopes = intersections.BCE_y(i,j-1) - positive_slope*...
                    intersections.BCE_x(i,j-1);
                reflected_characteristic = positive_slope * x_mesh+...
                                           diff_positive_slopes;
            end
        else
            negative_slope = -Negative_Slope(Wave.BCE_mu(i,j), Wave.BCE_phi(i,j));
            if j == i
                positive_slope = 0;
            else
                positive_slope = Positive_Slope(Wave.BCE_mu(i,j), Wave.BCE_phi(i,j));
            end
            incoming_characteristic  = intersections.BCE_y(i-1,j) + negative_slope*...
                                       (x_mesh - intersections.BCE_x(i-1,j));
            reflected_characteristic = intersections.BCE_y(i,j-1) + positive_slope*...
                                       (x_mesh - intersections.BCE_x(i,j-1));
        end
        y_i = find(diff(sign(incoming_characteristic - reflected_characteristic)));
        intersections.BCE_y(i,j) = reflected_characteristic(y_i(end));
        intersections.BCE_x(i,j) = x_mesh(y_i(end));
            
        x_i = find((x_mesh >= intersections.BCE_x(i,j)));
        if i == 1
            gammas(j,1:y_i(1)) = h + negative_slope * x_mesh(y_i(1));
            if j == 1
                prev_x_idx_reflected = find(x_mesh >= intersections.BCE_x(i,end));
            else
                prev_x_idx_reflected = find(x_mesh >= intersections.BCE_x(i,j-1));
            end
        elseif i == j
            prev_x_idx = find((x_mesh >= intersections.BCE_x(i-1,j)));
            gammas(j,prev_x_idx(1,1):end) = intersections.BCE_y(i-1,j) +...
                                            negative_slope *...
                                            (x_mesh(prev_x_idx(1,1:end)) -...
                                            x_mesh(prev_x_idx(1,1)));
        else
            prev_x_idx = find(x_mesh >= intersections.BCE_x(i-1,j));
            prev_x_idx_reflected = find(x_mesh >= intersections.BCE_x(i,j-1));
            gammas(j,prev_x_idx(1,1):end) = intersections.BCE_y(i-1,j) + ...
            negative_slope * (x_mesh(prev_x_idx(1,1):end) - x_mesh(prev_x_idx(1,1)));
            gammas(i,prev_x_idx_reflected(1,1):end) = ...
            intersections.BCE_y(i,j-1) + positive_slope *...
            (x_mesh(prev_x_idx_reflected(1,1):end) - x_mesh(prev_x_idx_reflected(1,1)));
        end
    end
end
for i = 1:N
    x_i = find(x_mesh <= intersections.BCE_x(1,end-(i-1)));
    negative_slope = (intersections.BCE_y(1,end-(i-1))-1)/...
                     (intersections.BCE_x(1,end-(i-1)))*x_range/x_steps;
    gammas(end-(i-1),1:length(x_i)) = Linear_Function(negative_slope,1,x_i,0);
    if i == N
    else
        x_i = find(x_mesh >= intersections.BCE_x(1,i) &...
                   x_mesh <= intersections.BCE_x(1,i+1));
        positive_slope = ((intersections.BCE_y(1,i+1)-intersections.BCE_y(1,i))...
        / (intersections.BCE_x(1,i+1)-intersections.BCE_x(1,i))) * x_range/x_steps;
        gammas(1, x_i(1):x_i(end)) = Linear_Function(positive_slope,...
                                     intersections.BCE_y(1,i), x_i - x_i(1), 0);
    end
end

%%%                           REGION DFG                                %%%
%%%=====================================================================%%%
jet_boundary = h + tan(phi_B) * x_mesh;
for i = 1:N
    for j = i:N
        diff_slopes_positive = 0;
        diff_slopes_negative = 0;
        if i == 1
            if j == 1
            diff_slopes = intersections.DFG_y(i,end) - Wave.BCE_alpha(j,end)...
                          *intersections.DFG_x(i,end);
            else
            diff_slopes = intersections.DFG_y(i,j-1) - Wave.BCE_alpha(j,end)...
                          *intersections.DFG_x(i,j-1);
            end
            incoming_characteristic = Linear_Function(Wave.BCE_alpha(j,end),...
            intersections.BCE_y(j,end), x_mesh, intersections.BCE_x(j,end));
            positive_slope = Positive_Slope(Wave.DFG_mu(i,j), Wave.DFG_phi(i,j));
            negative_slope = -Negative_Slope(Wave.DFG_mu(i,j), Wave.DFG_phi(i,j));
            if j == i
                reflected_characteristic = jet_boundary;
            else
                negative_slope = -Negative_Slope(Wave.DFG_mu(i,j), Wave.DFG_phi(i,j));
                diff_slopes_negative = intersections.DFG_y(i,j-1) -...
                                        negative_slope * intersections.DFG_x(i,j-1);
                reflected_characteristic = Linear_Function(negative_slope,...
                                           diff_slopes_negative, x_mesh, 0);
            end
        else
            incoming_characteristic = 0;
            reflected_characteristic = 0;
            negative_slope = -Negative_Slope(Wave.DFG_mu(i,j), Wave.DFG_phi(i,j));
            positive_slope = Positive_Slope(Wave.DFG_mu(i,j), Wave.DFG_phi(i,j));
            if j ~= i
                diff_slopes_positive = intersections.DFG_y(i-1,j)...
                - positive_slope * intersections.DFG_x(i-1,j);
                diff_slopes_negative = intersections.DFG_y(i,j-1)...
                - negative_slope * intersections.DFG_x(i,j-1);
                incoming_characteristic = Linear_Function(positive_slope,...
                                          diff_slopes_positive, x_mesh, 0);
                reflected_characteristic = Linear_Function(negative_slope,...
                                            diff_slopes_negative, x_mesh, 0);
            else
                diff_slopes_positive = intersections.DFG_y(i-1,j)...
                - positive_slope * intersections.DFG_x(i-1,j);
                incoming_characteristic = Linear_Function(positive_slope,...
                                          diff_slopes_positive, x_mesh, 0);
                reflected_characteristic = jet_boundary;
            end
        end
        y_i = find(diff(sign(incoming_characteristic - reflected_characteristic)));
        intersections.DFG_y(i,j) = incoming_characteristic(y_i(end));
        intersections.DFG_x(i,j) = x_mesh(y_i(end));
        
        x_idx = find(x_mesh >= intersections.DFG_x(i,j));
        if i == 1
            prev_Wave.x_idx = find(x_mesh >= intersections.BCE_x(j,end));
            if j == 1
                prev_x_idx_reflected = find(x_mesh >= intersections.DFG_x(i,end));
            else
                prev_x_idx_reflected = find(x_mesh >= intersections.DFG_x(i,j-1));
            end
            gammas(j,prev_Wave.x_idx(1,1):x_idx(1,1)) = intersections.BCE_y(j,end)...
            + positive_slope * (x_mesh(prev_Wave.x_idx(1,1):x_idx(1,1))...
            - x_mesh(prev_Wave.x_idx(1,1)));
            if j > 1
                gammas(j,prev_Wave.x_idx(1,1):x_idx(1,1)) = intersections.BCE_y(j,end)...
                + Wave.BCE_alpha(j,end) * (x_mesh(prev_Wave.x_idx(1,1):x_idx(1,1))...
                - x_mesh(prev_Wave.x_idx(1,1)));
                gammas(i,prev_x_idx_reflected(1,1):x_idx(1,1) + 1) =...
               intersections.DFG_y(i,j-1) + negative_slope * ...
               (x_mesh(prev_x_idx_reflected(1,1):x_idx(1,1)+1) -...
               x_mesh(prev_x_idx_reflected(1,1)));
            else
                diff_slopes = jet_boundary(x_idx(1,1)) - x_mesh(x_idx(1,1))...
                * tan(Wave.DFG_phi(i,j));
                jet_boundary(x_idx(1,1):end) = tan(Wave.DFG_phi(i,j))...
                * x_mesh(x_idx(1,1):end) + diff_slopes;
            end
        elseif j == i
            prev_x_idx = find(x_mesh >= intersections.DFG_x(i-1,j));
            gammas(j,prev_x_idx(1,1):x_idx(1,1) +1) = intersections.DFG_y(i-1,j)...
            + positive_slope * (x_mesh(prev_x_idx(1,1):x_idx(1,1) +1)...
            - x_mesh(prev_x_idx(1,1)));
            diff_slopes = jet_boundary(x_idx(1,1)) - x_mesh(x_idx(1,1))...
            * tan(Wave.DFG_phi(i,j));
            jet_boundary(x_idx(1,1):end) = tan(Wave.DFG_phi(i,j))...
            * x_mesh(x_idx(1,1):end) + diff_slopes;
        else
            prev_x_idx = find(x_mesh >= intersections.DFG_x(i-1,j));
            idx = find(x_mesh >= intersections.DFG_x(i,j));
            prev_x_idx_reflected = find(x_mesh >= intersections.DFG_x(i,j-1));
            gammas(i,prev_x_idx_reflected(1,1):end) = negative_slope *...
            (x_mesh(prev_x_idx_reflected(1,1):end)) + diff_slopes_negative; 
            gammas(j,prev_x_idx(1,1):end) = positive_slope *...
            (x_mesh(prev_x_idx(1,1):end)) + diff_slopes_positive;
        end
        if j == N - 1
            gammas(i,x_idx(1,1):end) = intersections.DFG_y(i,j)...
            + Wave.DFG_alpha(i,j) * (x_mesh(x_idx(1,1):end) - x_mesh(x_idx(1,1)));
        end
    end
end

%%%                           REGION HIJ                                %%%
%%%=====================================================================%%%
for i = 1:N
    for j = i:N
        diff_slopes_positive = 0;
        diff_slopes_negative = 0;
        positive_slope = 0;
        negative_slope = 0;
        if i == 1
            negative_slope = Wave.DFG_alpha(j,end);
            diff_slopes_negative = intersections.DFG_y(j,end) ...
            - negative_slope * intersections.DFG_x(j,end);
            incoming_characteristic = diff_slopes_negative + negative_slope...
            * x_mesh;
            if j == 1
                reflected_characteristic = 0 * x_mesh;
            else
                positive_slope = Positive_Slope(Wave.HIJ_mu(i,j-1),...
                                 Wave.HIJ_phi(i,j-1));
                diff_slopes_positive = intersections.HIJ_y(i,j-1)...
                - positive_slope * intersections.HIJ_x(i,j-1);
                reflected_characteristic = positive_slope...
                * x_mesh + diff_slopes_positive;
            end
        else
            negative_slope = -Negative_Slope(Wave.HIJ_mu(i-1,j), Wave.HIJ_phi(i-1,j));
            if j == i
                reflected_characteristic = 0 * x_mesh;
            else
                positive_slope = Positive_Slope(Wave.HIJ_mu(i,j), Wave.HIJ_phi(i,j));
                diff_slopes_positive = intersections.HIJ_y(i,j-1)...
                - positive_slope * intersections.HIJ_x(i,j-1);
                reflected_characteristic = intersections.HIJ_y(i,j-1)...
                + positive_slope * (x_mesh - intersections.HIJ_x(i,j-1));
            end
            incoming_characteristic = intersections.HIJ_y(i-1,j) +...
            negative_slope * (x_mesh - intersections.HIJ_x(i-1,j));
        end
        y_i = find(diff(sign(incoming_characteristic - reflected_characteristic)));
        intersections.HIJ_y(i,j) = reflected_characteristic(y_i);
        intersections.HIJ_x(i,j) = x_mesh(y_i); 
        
        x_idx = find(x_mesh >= intersections.HIJ_x(i,j));
        if i == 1
            prev_Wave.x_idx = find(x_mesh >= intersections.DFG_x(j,end));
            gammas(j,prev_Wave.x_idx(1,1):end) = diff_slopes_negative...
            + negative_slope * x_mesh(prev_Wave.x_idx(1,1):end);
            if j > 1
                prev_x_idx_reflected = find(x_mesh >= intersections.HIJ_x(i,j-1));
                gammas(i,prev_x_idx_reflected(1,1):end) =...
                intersections.HIJ_y(i,j-1) + positive_slope...
                *(x_mesh(prev_x_idx_reflected(1,1):end)-...
                x_mesh(prev_x_idx_reflected(1,1)));
            end
        elseif j == i
            prev_x_idx = find(x_mesh >= intersections.HIJ_x(i-1,j));
            gammas(j,prev_x_idx(1,1):end) = intersections.HIJ_y(i-1,j)...
            + negative_slope * (x_mesh(prev_x_idx(1,1):end)...
            - x_mesh(prev_x_idx(1,1)));
        else
            prev_x_idx = find(x_mesh >= intersections.HIJ_x(i-1,j));
            prev_x_idx_reflected = find(x_mesh >= intersections.HIJ_x(i,j-1));
            idx = find(x_mesh >= intersections.HIJ_x(i,j));
            gammas(j,prev_x_idx(1,1):end) = intersections.HIJ_y(i-1,j)...
            + negative_slope * (x_mesh(prev_x_idx(1,1):end)-...
            x_mesh(prev_x_idx(1,1)));
            gammas(i,prev_x_idx_reflected(1,1):end) = intersections.HIJ_y(i,j-1)...
            + positive_slope * (x_mesh(prev_x_idx_reflected(1,1):end)...
            - x_mesh(prev_x_idx_reflected(1,1)));
            if j == N-1
                gammas(i,x_idx(1,1):end) = intersections.HIJ_y(i,j)...
                + Wave.HIJ_alpha(i,j) * (x_mesh(x_idx(1,1):end)...
                - x_mesh(x_idx(1,1)));
            end
        end
    end
end

%%%                      MINOR IRRELEVANT FIX                           %%%
x_i = find(x_mesh >= intersections.HIJ_x(N,N));
gammas(N,x_i(1):end) = -gammas(N,x_i(1):end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       Creating heat map                             %%%
%%%     function: Create_Heat_Map_Geometry constructs all patches       %%%
%%%              where angels of attack will be plotted                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap = Create_Heat_Map_Geometry(intersections, N, x_range, jet_boundary);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Creating Streamline for given initial height              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_height_streamline_1 = 0.5;
initial_height_streamline_2 = 0.75;
initial_height_streamline_3 = 0.25;

[streamline_1, press_dist_1, x_crossings_1] =...
Create_Streamline(initial_height_streamline_1, x_mesh, mach_initial,...
gamma, pressure_total, N, x_steps, intersections, gammas, Wave);
[streamline_2, press_dist_2, x_crossings_2] =...
Create_Streamline(initial_height_streamline_2, x_mesh, mach_initial,...
gamma, pressure_total, N, x_steps, intersections, gammas, Wave);
[streamline_3, press_dist_3, x_crossings_3] =...
Create_Streamline(initial_height_streamline_3, x_mesh, mach_initial,...
gamma, pressure_total, N, x_steps, intersections, gammas, Wave);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       Find shock location                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_H_2 = find(x_mesh > intersections.HIJ_x(2,2));
idx_H_3 = find(x_mesh > intersections.HIJ_x(3,3));
% gamma_H_2 = gammas(2,idx_H_2(1):end);
% gamma_H_3 = gammas(3,idx_H_3(1):end);

[shock_x,shock_y] = polyxpoly(x_mesh(idx_H_2:end),gammas(2,idx_H_2:end),x_mesh(idx_H_3:end),gammas(3,idx_H_3:end));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        Plotting section                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boundary.side_x = [0 , 0]; boundary.side_y = [0, 1]; 
boundary.end_x = [0, x_range]; boundary.end_y = [0, 0];

% %%%                       MACH DISTRIBUTION                             %%%
% %%%=====================================================================%%%
% figure(1)
% set(gcf,'Position',[200 200 1500 600])
% hold on; box on;
% xlabel('x [m]'); ylabel('y [m]');
% xlim([0 intersections.HIJ_x(2,end)])
% ylim([0 max(jet_boundary)+0.1])
% plot(x_mesh, gammas, 'b')
% plot(boundary.side_x, boundary.side_y, 'k')
% plot(boundary.end_x, boundary.end_y, 'k')
% plot(x_mesh,jet_boundary,'k')
% %%                        COLORMAPPING                                 %%%
% patch(colormap.AB0_x,colormap.AB0_y,mach_initial,'EdgeColor','none')
% for i = 1:N-1
%     patch(colormap.ABC_x{i},colormap.ABC_y{i},Wave.ABC_mach(i),'EdgeColor','none')
% end
% patch(colormap.ACD_x,colormap.ACD_y,mach_B,'EdgeColor','none')
% for i = 1:N-1
%     for j = i:N-1
%         patch(colormap.BCE_x{i,j},colormap.BCE_y{i,j},Wave.BCE_mach(i,j),'EdgeColor','none')
%     end
% end
% for i = 1:N-1
%     patch(colormap.CDEF_x{i}, colormap.CDEF_y{i}, Wave.BCE_mach(i,end),'EdgeColor','none')
% end
% patch(colormap.EFH_x, colormap.EFH_y, Wave.BCE_mach(end,end),'EdgeColor','none')
% for i = 1:N-1
%     for j = i:N-1
%         patch(colormap.DFG_x{i,j}, colormap.DFG_y{i,j}, Wave.DFG_mach(i,j),'Edgecolor','none')
%     end
% end
% for i = 1:N-1
%     patch(colormap.GFHI_x{i}, colormap.GFHI_y{i}, Wave.DFG_mach(i,end),'EdgeColor','none')
% end
% for i = 1:N-1
%     for j = i:N-1
%         patch(colormap.HIJ_x{i,j},colormap.HIJ_y{i,j},Wave.HIJ_mach(i,j),'EdgeColor','none')
%     end
% end
% patch(colormap.IJK_x, colormap.IJK_y, Wave.DFG_mach(end,end),'Edgecolor','none')
% c = colorbar('Limits',[2,3]);
% colormap;

%%%                           STREAMLINES                               %%%
%%%=====================================================================%%%
% figure(2)
% set(gcf,'Position',[100 100 1500 600])
% hold on; box on;
% xlabel('x [m]'); ylabel('y [m]');
% xlim([0 x_range])
% ylim([0 max(jet_boundary)+0.1])
% plot(x_mesh, gammas, 'b')
% plot(boundary.side_x, boundary.side_y, 'k')
% plot(boundary.end_x, boundary.end_y, 'k')
% plot(x_mesh,jet_boundary,'k')
% % plot(x_mesh,streamline_1,'k')
% plot(x_mesh,streamline_2,'k')
% % plot(x_mesh,streamline_3,'k')


%%%               PRESSURE DISTRIBUTION ALONG STREAMLINE                %%%
%%%=====================================================================%%%
% figure(3)
% set(gcf,'Position',[300 300 1000 450])
% box on; grid on;hold on;
% xlim([0 13.3])
% xlabel('x [m]'); ylabel('p/p_{tot} [-]');
% plot(x_crossings_1,  press_dist_1,'DisplayName','Streamline 1 at H/2','LineStyle','-','Color','k')
% plot(x_crossings_2,  press_dist_2,'DisplayName','Streamline 2 at 3H/4','LineStyle','--','Color','k')
% legend

%%%                         EXPANSION FAN                               %%%
%%%=====================================================================%%%
% expansion_fan = zeros(length(alpha_list),length(x_mesh));
% for i = 1:N
% expansion_fan(i,:) = Linear_Function(alpha_list(i),1,x_mesh,0);
% end
% figure(5)
% set(gcf,'Position',[200 200 800 400])
% hold on; box on;
% xlabel('x [m]'); ylabel('y [m]');
% ylim([0 1])
% xlim([0 5])
% plot(x_mesh,expansion_fan(1,:),'k');
% plot(x_mesh,expansion_fan(2:end-1,:),'k--');
% plot(x_mesh,expansion_fan(end,:),'k');
% % legend('Characteristic')

%%%                         CHARACTERISTICS                             %%%
%%%=====================================================================%%%
figure(6)
set(gcf,'Position',[100 100 1500 600])
hold on; box on
xlabel('x [m]'); ylabel('y [m]');
xlim([0 x_range])
ylim([0 max(jet_boundary)+0.1])
plot(x_mesh, gammas, 'b')
plot(boundary.side_x, boundary.side_y, 'k')
plot(boundary.end_x, boundary.end_y, 'k')
plot(shock_x,shock_y,'ro')
plot(x_mesh,jet_boundary,'k')

%%                        SHOCK FORMATION                             %%%
%%%=====================================================================%%%
% Mach_list  = [1.70 1.75 1.80 1.85 1.90 1.95 2.00 2.05 2.10];
% shock_list = [9.7758 9.9051 10.3506 10.7506 11.0227 11.3414 11.7069 11.7261 12.2392];
% figure(7)
% set(gcf,'Position',[200 200 800 400])
% ylim([1.65 2.15])
% hold on; box on; grid on;
% ylabel('Exit Mach'); xlabel('Shock Location [m]');
% plot(shock_list, Mach_list,'k-o')







