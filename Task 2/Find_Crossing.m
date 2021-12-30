function [lowest_idx, gamma] = Find_Crossing(streamline_indices, prev_idx, x_steps, N)
%FIND_CROSSING Summary of this function goes here
%   Detailed explanation goes here
lowest_idx = x_steps;
gamma = 0;
for i = 1:N
    if  (streamline_indices.index(i) > prev_idx && streamline_indices.index(i) < lowest_idx)
        lowest_idx = streamline_indices.index(i);
        gamma = streamline_indices.gamma(i);
    end
end

