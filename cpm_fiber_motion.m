function [out, tt, traj_x, traj_y] = cpm_fiber_motion(sigma, fiber,J_cc, J_cm, J_cf, J_p, alpha, temp_po, lam_area, cell_size, TEMP, max_ite)
% an updated way of running the MC simulation. The main difference is the
% random sampling is from the boundary points which is recorded and later
% updated. This algorithm is about 6 times fater than the original one
% J_p is the strength of motion energy, alpha is the probability of keeping
% the current polarity. (1-alpha) is the probability of choosing a new dir.

Nx = size(sigma, 1);
Ny = size(sigma, 2);

tt = [];
traj_x = [];
traj_y = [];

max_step = Nx * Ny * max_ite;

% the current cell area
temp_area = sum(sum(sigma));

[c_bound_x,c_bound_y] = detect_boundary(sigma, 1); % record the boundary of cells in cell array
[m_bound_x, m_bound_y] = detect_boundary(sigma, 0); % record the boundary of medium in cell array

[xx,yy] = meshgrid(1:size(sigma,1),1:size(sigma,2));


scc_try = 0;
f_record = 1;
% rand_theta = rand * 2 * pi;
rand_theta = 3 *pi / 4;
p_0 = [cos(rand_theta), sin(rand_theta)];

for step = 1 : max_step
    % get all the possible sample points at boundaries
    sample_bound_x = [c_bound_x, m_bound_x];
    sample_bound_y = [c_bound_y, m_bound_y];
    
    % find the point we want to flip
    rand_idx = randi(length(sample_bound_x));
    source_x = sample_bound_x(rand_idx);
    source_y = sample_bound_y(rand_idx);
    source = sigma(source_x, source_y);
    f_ind = fiber(source_x, source_y);
    
    % find all the neibors and the target we want the source pst flip to
    neibor_x = [mod(source_x-2, Nx)+1, source_x, mod(source_x, Nx)+1];
    neibor_y = [mod(source_y-2, Ny)+1, source_y, mod(source_y, Ny)+1];
    neibors = sigma(neibor_x, neibor_y); neibors = reshape(neibors,[9,1]); 
    neibors(5) = [];
    target = neibors(randi(8)); 

    if source ~= target
        scc_try = scc_try + 1;
        % energy is calculated as after - before flipping 
        % interfacial energy is calculated for all the neibors
        nb_c = sum(neibors == 1); nb_m = sum(neibors == 0);
        if source == 1
            sfc_energy = J_cm * nb_c - (J_cc * nb_c + J_cm * nb_m);
        else
            sfc_energy = J_cc * nb_c + J_cm * nb_m - J_cm * nb_c;
        end
        
        % need periodic boundary condition
        tot_mass = sum(sigma(:));
        
%         xc_sin = sigma / tot_mass .*sin(xx*pi/Nx);
%         xc_cos = sigma / tot_mass .*cos(xx*pi/Nx);
%         mc_x = atan2(sum(xc_sin(:)), sum(xc_cos(:))) / pi * Nx;
%         
%         yc_sin = sigma / tot_mass .*sin(yy*pi/Ny);
%         yc_cos = sigma / tot_mass .*cos(yy*pi/Ny);
%         mc_y = atan2(sum(yc_sin(:)), sum(yc_cos(:))) / pi * Ny;
        mc = xx.*sigma;
        mc_x = sum(mc(:))/tot_mass;
        mc = yy.*sigma;
        mc_y = sum(mc(:))/tot_mass;
        ran1 = rand;
        if ran1 > alpha && mod(scc_try, 100) == 1
            p_0 = repolarize(sigma, mc_x, mc_y, temp_po);
        end
        p_energy = 0;
        
%         [vec_x, vec_y] = vec_dis(source_x, source_y, mc_x, mc_y, Nx, Ny);
        % notice here the energy function is different from the polarity
        % moduel which I don't understand why
        if source == 1 %lose energy
            p_energy = -J_p*((source_x-mc_x)*p_0(2) + (source_y-mc_y)*p_0(1));
        elseif source == 0 % gain energy
            p_energy = J_p*((source_x-mc_x)*p_0(2) + (source_y-mc_y)*p_0(1));
        end
            
%         sfc_energy = -J_cm * sum(neibors ~= source) + J_cm * sum(neibors ~= target);
        
        % cell-fiber energy is a potential energy, which only depends on
        % psts and should be easy to calculate
        cf_energy = 0;
        if source == 1
            area_energy = lam_area * (-(temp_area - cell_size)^2 + (temp_area - 1 - cell_size)^2);
            if f_ind == 1
                cf_energy = -J_cf;
            end
        elseif source == 0
            area_energy = lam_area * (-(temp_area - cell_size)^2 + (temp_area + 1 - cell_size)^2);
            if f_ind == 1
                cf_energy = J_cf;
            end
        end
        
        tot_energy = sfc_energy + area_energy + cf_energy + p_energy;
        if tot_energy > 0
            prob = exp(-tot_energy / TEMP);
        elseif tot_energy <= 0
            prob = 1;
        end
        
        if prob > rand
            % update area
            if source == 1
                temp_area = temp_area - 1;
            elseif source == 0
                temp_area = temp_area + 1;
            end
            % flip the source sigma to target
            sigma(source_x,source_y) = target;
            % update the boundary points
            for i = 1 : 3
                for j = 1 : 3
                    % cancel all the neiboring points
                    del_idx = (c_bound_x == neibor_x(i) & c_bound_y == neibor_y(j));
                    c_bound_x(del_idx) = [];
                    c_bound_y(del_idx) = [];
                    del_idx = (m_bound_x == neibor_x(i) & m_bound_y == neibor_y(j));
                    m_bound_x(del_idx) = [];
                    m_bound_y(del_idx) = [];
                    if isboundary(neibor_x(i),neibor_y(j),sigma)
                        % join the newly formed neiboring points
                        if sigma(neibor_x(i), neibor_y(j)) == 1
                            c_bound_x(end+1) = neibor_x(i);
                            c_bound_y(end+1) = neibor_y(j);
                        else
                            m_bound_x(end+1) = neibor_x(i);
                            m_bound_y(end+1) = neibor_y(j);
                        end
                    end
                end
            end
        end
        % shift to center
        if max(c_bound_x) > Nx -5 || min(c_bound_x) < 5 || max(c_bound_y) > Ny -5 || min(c_bound_y) < 5
            sigma = circshift(sigma,[floor(Ny/2-mc_y),floor(Nx/2-mc_x)]);
            [c_bound_x,c_bound_y] = detect_boundary(sigma, 1); 
            [m_bound_x, m_bound_y] = detect_boundary(sigma, 0);
        end
       
    end
    
    if mod(scc_try, 100) == 1
        tt(end+1) = scc_try;
        traj_x(end+1) = mc_x;
        traj_y(end+1) = mc_y;
    end
    
    % record movie
    if mod(step, Nx*Ny) == 0
        clf;
        imagesc(sigma); axis square; set(gca,'YDir','normal')
        hold on
        [fx,fy]=find(fiber);
        plot(fy,fx,'ro');
        plot([mc_x,mc_x+5*p_0(1)],[mc_y,mc_y+5*p_0(2)],'k','linewidth',2);
        im(f_record) = getframe(gcf);
        f_record = f_record + 1;
    end
end
% v = VideoWriter('time_evolve'); v.FrameRate = 5;
% open(v); writeVideo(v,im); close(v);

out = sigma;
% scc_try