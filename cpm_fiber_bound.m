function out = cpm_fiber_bound(sigma, fiber,J_cc, J_cm, J_cf, lam_area, cell_size, TEMP, max_ite)
% an updated way of running the MC simulation. The main difference is the
% random sampling is from the boundary points which is recorded and later
% updated. This algorithm is about 6 times fater than the original one

Nx = size(sigma, 1);
Ny = size(sigma, 2);

max_step = Nx * Ny * max_ite;

% the current cell area
temp_area = sum(sum(sigma));

[c_bound_x,c_bound_y] = detect_boundary(sigma, 1); % record the boundary of cells in cell array
[m_bound_x, m_bound_y] = detect_boundary(sigma, 0); % record the boundary of medium in cell array



scc_try = 0;

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
        
        tot_energy = sfc_energy + area_energy + cf_energy;
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
        
    end
end

out = sigma;
% scc_try