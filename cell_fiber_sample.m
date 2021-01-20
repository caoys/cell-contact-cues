% the sample program using latin hyperbolic sampling on J_cc, J_cf, and
% dphi
close all
clear

Nx = 50;
Ny = 50;
J_cc = 0;
% J_cm = 10; % cell-medium adhesion
%J_cf = -50; % cell-fiber adhesion
lam_area = 10; % cell-size constraint
cell_size = 15^2; % cell size
TEMP = 10;
N_fiber = 10; % the total fiber numbers in simulation box which is determined by the fiber density and cell size

% max MC steps in Nx*Ny box
max_it = 100;

% get the max number of cells
max_trial = 2000;
% fiber_var = 60;
%J_cf = 30;
tic
sample_no = 1;
cell_fiber_orient = [];
cell_aspect = [];
cell_contract = []; % record J_cc
cell_fiber_var = []; % record fiber_var
cell_adh = []; % record J_cf
cell_seed =[];

ran_p = lhsdesign(max_trial, 4);

sigma_init = zeros(Nx,Ny);
cell_length = sqrt(cell_size);
sigma_init(Nx / 2 - floor(cell_length / 2) : Nx / 2 - floor(cell_length / 2) + cell_length-1, Ny / 2 - floor(cell_length / 2) : Ny / 2 - floor(cell_length / 2) + cell_length-1) = 1;
temp_area = sum(sum(sigma_init));

% fiber step length
fl = 5;
sigma_l = 2;
density = 900;

gamma = 3;
% rg = 10;
while sample_no <= max_trial
% for J_cf = 10: 10 : 20
%for fiber_var = 70 : 10 : 80


    % sample J_cf~(10,60), fiber_var~(10,80)
    J_cm = ran_p(sample_no,1) * 20 + 10;
    J_cf = ran_p(sample_no,2) * gamma * J_cm;
    seed_N = ceil(ran_p(sample_no, 3) * 17); % fiber seed numbers
    
    fiber_var = 2.5 + (seed_N - 1) * 5 + ran_p(sample_no,4)*5; % fiber variance
    rn = round(density / fl / (seed_N + 1)); % steps
    
    rand_angle = 0;%(2 * rand - 1) * pi / 2;
    [fiber, dir] = random_walk_fiber(seed_N, rn, fl, sigma_l, fiber_var/180*pi, 0, Nx, Ny);
    
%     fiber_var = ran_p(sample_no,3) * 88 + 1;

    % the orientation difference/variance of fibers
%     fiber_var = 90;
%     tic
%     for k = 1 : max_trial
        % choose a random mean orietation 
        sigma = sigma_init;
%         [fiber, dir] = random_fiber(rand_angle, fiber_var/180*pi, rg, Nx, Ny, N_fiber);
        sigma = cpm_fiber_bound(sigma,fiber, -J_cc,J_cm,-J_cf,lam_area,cell_size,TEMP,max_it);
        [cm_x, cm_y, ecc, theta, maj_ax, min_ax] = cell_analysis(sigma);
        % the x-axis for theta is on the negative direction
        cell_fiber_orient(end+1) = (theta-dir)/180*pi - rand_angle;
        cell_seed(end+1) = seed_N;
        cell_aspect(end+1) = maj_ax / min_ax;
        cell_contract(end+1) = J_cm;
        cell_fiber_var(end+1) = fiber_var;
        cell_adh(end+1) = J_cf;
%     end
%     toc
% end
    sample_no = sample_no + 1;
end
save(strcat('new_sample_gamma_',int2str(gamma),'_den_',int2str(density),'.mat'),'cell_fiber_var','cell_adh','cell_contract','cell_seed','cell_fiber_orient','cell_aspect');
toc
