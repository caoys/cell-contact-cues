function [xl,yl] = line_int(x0,y0,xend,yend)
% find all the integer grid points determined by starting point of (x0,y0)
% and ending point (xend, yend). All poits are stored in xl and yl,
% respectively

k = (yend - y0) / (xend - x0);
if isnan(k)
    k = 100;
end

kinv = (xend - x0) / (yend - y0);
if isnan(kinv)
    kinv = 100;
end

xg = floor(x0) : 1 : floor(xend);
yg = floor(y0) : 1 : floor(yend);

if k <=1 
    xl = xg;
    yl = floor(k*(xl-x0) + y0);
else
    yl = yg;
    xl = floor(kinv *(yl-y0) + x0);
end


