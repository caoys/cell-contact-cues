function [out, dir] = random_walk_fiber(N, rn, fl, sigma_l, delta, phi0, Nx, Ny)
% algorithm to generate fibers based on random walk growth with spatial
% length of fl, drawn from normal distribution with sigma_l and each step
% with angle sampled from uniform distribution
% of delta, and N+1 seeds. One seed is always located at the center of
% (Nx/2,Ny/2)
% rn*fl*(N+1) should be the total fiber occupancy

out = zeros(Nx, Ny);
dir = phi0;
% generate all spatial points
xs = zeros(rn, N+1);
ys = zeros(rn, N+1);

xs(1,1) = Nx / 2;
ys(1,1) = Ny / 2;

% initial starting point
xs(1, 2 : end) = randi([1, Nx], [1,N]);
ys(1, 2 : end) = randi([1, Ny], [1,N]);

xg = [];
yg = [];
% repeat/grow for rn times
for i = 1 : rn - 1
    % get random directions and random length
    step_l = abs(normrnd(fl, sigma_l,[1, N+1]));
    step_phi = unifrnd(-delta+phi0, delta+phi0,[1, N+1]);
    dir = dir + sum(step_phi);
    xs(i+1,:) = xs(i, :) + step_l./sqrt(1+(tan(step_phi)).^2);
    ys(i+1,:) = ys(i, :) + tan(step_phi).*step_l./sqrt(1+(tan(step_phi)).^2);
    for j = 1 : N+1
        [xn, yn] = line_int(xs(i,j),ys(i,j), xs(i+1,j),ys(i+1,j));
        xg = [xg,xn];
        yg = [yg,yn];
    end
end

m_idx = sub2ind(size(out),mod(xg, Nx)+1,mod(yg,Ny)+1);
out(m_idx) = 1;
out = out';
dir = dir / rn / (N+1);


