function p_out = repolarize(sigma, mc_x, mc_y, temp_po)
% a function reselect a new polarity based on 
% P(p_0) = int (p*p_0)
p_theta = linspace(0, 2*pi, 20);
TEMP = temp_po;
% tot_mass = sum(sigma(:));
% [xx,yy] = ndgrid(1:size(sigma,1),1:size(sigma,2));
% mc_x = sum(xx(:).*sigma(:))/tot_mass;
% mc_y = sum(yy(:).*sigma(:))/tot_mass;

[xx,yy] = meshgrid(1:size(sigma,1),1:size(sigma,2));

% sigma_s = circshift(sigma, [floor(Ny/2-mc_y),floor(Nx/2-mc_x)]);

prob = zeros(size(p_theta));
for k = 1 : length(p_theta)
%     prob(k) = 0;
%     for i = 1 : Nx
%         for j = 1 : Ny
%             if sigma(i,j) == 1
%                 [vec_x, vec_y] = vec_dis(i, j, mc_x, mc_y, Nx, Ny);
%                 prob(k) = prob(k) + abs(vec_x*cos(p_theta(k)) + vec_y*sin(p_theta(k)));
%             end
%         end
%     end
    %get the intenral energy based on p
    u_internal = -abs(sigma.*(xx-mc_x)*cos(p_theta(k)) + sigma.*(yy-mc_y)*sin(p_theta(k)));
    prob(k) = exp(-sum(u_internal(:))/TEMP);
end
prob = prob / sum(prob(:));
atemp = 0;
b = rand;
for q = 1 : length(prob)
    if b > atemp && b < atemp+prob(q)
        miu = q;
        break;
    end
    atemp = atemp + prob(q);
end
p_out = [cos(p_theta(miu)), sin(p_theta(miu))];