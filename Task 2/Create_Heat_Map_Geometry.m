function [colormap] = Create_Heat_Map_Geometry(intersections, N, x_range, jet_boundary)
%CREATE_HEAT_MAP_GEOMETRY Creates triangles and squares where local mach
%values can be plotted

%%%                          REGION OUTLET                              %%%
%%%=====================================================================%%%
colormap.AB0_x = [0; 0; intersections.BCE_x(1,1)];
colormap.AB0_y = [0; 1; intersections.BCE_y(1,1)];

%%%                           REGION ABC                                %%%
%%%=====================================================================%%%
for i = 1:N-1
colormap.ABC_x{i} = [0; intersections.BCE_x(1,i); intersections.BCE_x(1,i+1)];
colormap.ABC_y{i} = [1; intersections.BCE_y(1,i); intersections.BCE_y(1,i+1)];
end

%%%                           REGION ACD                                %%%
%%%=====================================================================%%%
colormap.ACD_x = [0; intersections.BCE_x(1,end); intersections.DFG_x(1,1)];
colormap.ACD_y = [1; intersections.BCE_y(1,end); intersections.DFG_y(1,1)];

%%%                           REGION BCE                                %%%
%%%=====================================================================%%%
for i = 1:N-1
    for j = i:N-1
        if j == i
        colormap.BCE_x{i,j} = [intersections.BCE_x(i,j);...
        intersections.BCE_x(i,j+1); intersections.BCE_x(i+1,j+1)];
        colormap.BCE_y{i,j} = [intersections.BCE_y(i,j);...
        intersections.BCE_y(i,j+1); intersections.BCE_y(i+1,j+1)];
        else
        colormap.BCE_x{i,j} = [intersections.BCE_x(i,j);...
        intersections.BCE_x(i,j+1); intersections.BCE_x(i+1,j+1);...
        intersections.BCE_x(i+1,j)];
        colormap.BCE_y{i,j} = [intersections.BCE_y(i,j);...
        intersections.BCE_y(i,j+1); intersections.BCE_y(i+1,j+1);...
        intersections.BCE_y(i+1,j)];
        end
    end
end

%%%                           REGION CDEF                               %%%
%%%=====================================================================%%%
for i = 1:N-1
    colormap.CDEF_x{i} = [intersections.BCE_x(i,end); intersections.BCE_x(i+1,end);...
                          intersections.DFG_x(1,i+1); intersections.DFG_x(1,i)];
    colormap.CDEF_y{i} = [intersections.BCE_y(i,end); intersections.BCE_y(i+1,end);...
                          intersections.DFG_y(1,i+1); intersections.DFG_y(1,i)];
end

%%%                           REGION EFH                                %%%
%%%=====================================================================%%%
colormap.EFH_x = [intersections.BCE_x(end,end); intersections.DFG_x(1,end);...
                  intersections.HIJ_x(1,1)];
colormap.EFH_y = [intersections.BCE_y(end,end); intersections.DFG_y(1,end);...
                  intersections.HIJ_y(1,1)];

%%%                           REGION DFG                                %%%
%%%=====================================================================%%%
for i = 1:N-1
    for j = i:N-1
        if j == i
        colormap.DFG_x{i,j} = [intersections.DFG_x(i,j);...
        intersections.DFG_x(i,j+1); intersections.DFG_x(i+1,j+1)];
        colormap.DFG_y{i,j} = [intersections.DFG_y(i,j);...
        intersections.DFG_y(i,j+1); intersections.DFG_y(i+1,j+1)];
        else
        colormap.DFG_x{i,j} = [intersections.DFG_x(i,j);...
        intersections.DFG_x(i,j+1); intersections.DFG_x(i+1,j+1);...
        intersections.DFG_x(i+1,j)];
        colormap.DFG_y{i,j} = [intersections.DFG_y(i,j);...
        intersections.DFG_y(i,j+1); intersections.DFG_y(i+1,j+1);...
        intersections.DFG_y(i+1,j)];
        end
    end
end

%%%                           REGION GFHI                               %%%
%%%=====================================================================%%%
for i = 1:N-1
    colormap.GFHI_x{i} = [intersections.DFG_x(i,end);...
    intersections.DFG_x(i+1,end); intersections.HIJ_x(1,i+1);...
    intersections.HIJ_x(1,i)];
    colormap.GFHI_y{i} = [intersections.DFG_y(i,end);...
    intersections.DFG_y(i+1,end); intersections.HIJ_y(1,i+1);...
    intersections.HIJ_y(1,i)];
end

%%%                           REGION HIJ                                %%%
%%%=====================================================================%%%
for i = 1:N-1
    for j = i:N-1
        if j == i
        colormap.HIJ_x{i,j} = [intersections.HIJ_x(i,j);...
        intersections.HIJ_x(i,j+1); intersections.HIJ_x(i+1,j+1)];
        colormap.HIJ_y{i,j} = [intersections.HIJ_y(i,j);...
        intersections.HIJ_y(i,j+1); intersections.HIJ_y(i+1,j+1)];
        else
        colormap.HIJ_x{i,j} = [intersections.HIJ_x(i,j);...
        intersections.HIJ_x(i,j+1); intersections.HIJ_x(i+1,j+1);...
        intersections.HIJ_x(i+1,j)];
        colormap.HIJ_y{i,j} = [intersections.HIJ_y(i,j);...
        intersections.HIJ_y(i,j+1); intersections.HIJ_y(i+1,j+1);...
        intersections.HIJ_y(i+1,j)];
        end
    end
end

%%%                           REGION IJK                                %%%
%%%=====================================================================%%%
colormap.IJK_x = [intersections.HIJ_x(2,end); intersections.DFG_x(end,end);...
                  x_range];
colormap.IJK_y = [intersections.HIJ_y(2,end); intersections.DFG_y(end,end);...
                  jet_boundary(end)];
end

