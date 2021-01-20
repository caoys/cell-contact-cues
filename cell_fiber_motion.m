% program to run cpm motion
close all
clear

Nx = 100;
Ny = 100;
J_cm = 10; % cell-medium adhesion

lam_area = 10; % cell-size constraint

cell_size = 20^2; % cell size
TEMP = 10;
max_it = 100; % the true step is Nx*Ny*max_it

sigma = zeros(Nx,Ny);
cell_length = sqrt(cell_size);
sigma(Nx / 2 - floor(cell_length / 2) : Nx / 2 - floor(cell_length / 2) + cell_length-1, Ny / 2 - floor(cell_length / 2) : Ny / 2 - floor(cell_length / 2) + cell_length-1) = 1;
temp_area = sum(sum(sigma));

% N_fiber = 17;
fl = 5;
sigma_l = 2;
density = 900;

seed_N = 3; % from 1 to 17
fiber_var = 2.5 + (seed_N - 1) * 5;
rn = round(density / fl / (seed_N + 1));

tic

    J_cc = 10;
    J_cf = 20; % cell-fiber adhesion
    J_p = -10; % direction potential
    alpha = 0.5; % the persistence of polarity
    temp_po = 10^5;
%     fiber_var = 80;
%     fiber = random_fiber(0, pi*fiber_var/180, Nx, Ny, N_fiber);
    [fiber, dir] = random_walk_fiber(seed_N, rn, fl, sigma_l, fiber_var/180*pi, 0, Nx, Ny);
    [sigma, tt, xt, yt] = cpm_fiber_motion(sigma,fiber, -J_cc,J_cm, -J_cf, J_p, alpha, temp_po, lam_area,cell_size,TEMP,max_it);

    [cm_x, cm_y, ecc, theta, maj_ax, min_ax] = cell_analysis(sigma);
    % ecc
    orient = theta / 180 * pi;
    asp = maj_ax / min_ax;

    % save('traj','tt','xt','yt');

    % figure(1); imagesc(sigma); axis square; set(gca,'YDir','normal')
    % hold on
    % [fx,fy]=find(fiber);
    % plot(fy,fx,'ro');
    % savefig('shape.fig');
    save('cf_motion_3.mat','orient','asp','tt','xt','yt');

toc
