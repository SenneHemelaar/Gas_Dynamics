function [streamline_indices] = Streamline_Crossing(streamline, prev_idx, N, gammas)
%STREAMLINE_CROSSING
streamline_indices = struct('index',zeros(N,1),'gamma',zeros(N,1));
for i = 1:N
    crossing_index = find(diff(sign(streamline(prev_idx:end) - gammas(i,prev_idx:end))));
    if (isempty(streamline_indices) && size(crossing_index) > 0)
        streamline_indices.index(i) = prev_idx + min(crossing_index);
        streamline_indices.gamma(i) = i;
    elseif (size(crossing_index) > 0)
        streamline_indices.index(i) = prev_idx + min(crossing_index);
        streamline_indices.gamma(i) = i;
    end
end
end

