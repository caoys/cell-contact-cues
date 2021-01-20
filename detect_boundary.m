function [idx_x,idx_y] = detect_boundary(sigma, ele)
% this is the function of detecting boundaries of ele in matrix sigma.
% Boundary is defined as the for any (i,j) in sigma, if any neighor is
% different, (i,j) is a boundary point, and stored in (idx_x, idx_y)

Nx = size(sigma, 1);
Ny = size(sigma, 2);
idx_x = [];
idx_y = [];

for i = 1 : Nx
    for j = 1 : Ny
        if sigma(i,j) == ele
            if isboundary(i, j, sigma)
                idx_x(end + 1) = i;
                idx_y(end + 1) = j;
            end
        end
    end
end

% idx_x = cell2mat(idx_x);
% idx_y = cell2mat(idx_y);