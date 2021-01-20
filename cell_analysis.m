function [cm_x, cm_y, ecc, theta, maj_ax, min_ax]=cell_analysis(sigma)

Nx = size(sigma, 1);
Ny = size(sigma, 2);
[xx,yy] = meshgrid(1:size(sigma,1),1:size(sigma,2));
tot_mass = sum(sigma(:));
mc = xx.*sigma;
mc_x = sum(mc(:))/tot_mass;
mc = yy.*sigma;
mc_y = sum(mc(:))/tot_mass;
sigma = circshift(sigma,[floor(Ny/2-mc_y),floor(Nx/2-mc_x)]);

% cc = bwconncomp(sigma);
s = regionprops(sigma, 'Area','Centroid','Orientation','MajorAxisLength', ...
    'MinorAxisLength', 'Eccentricity'); 
idx = find([s.Area] == max([s.Area])); 

cm_x = s(idx).Centroid(1);
cm_y = s(idx).Centroid(2);
ecc = s(idx).Eccentricity;
theta = s(idx).Orientation;
maj_ax = s(idx).MajorAxisLength;
min_ax = s(idx).MinorAxisLength;

% Outlines=bwboundaries(cc);
% Xout=Outlines(:,1);
% Yout=Outlines(:,2);
% imagesc(cc);axis equal
% axis off
% hold on
% plot(Xout,Yout,'r');