function shooting
%SHOOTING (example) Convert BVP to IVP - guess init. condition & check if boundary conditions are met
%   Uses bisection to search for best init value for one other unknown variable

z1 = -3; z2 = 3; % bounds to search for unknown init value
yb = 2; % actual value of known variable at BC2
y2 = init2bound(z2);
num_iter = 1e6;
for i=1:num_iter
    znew = (z1+z2)/2;
    yb_new = init2bound(znew);
    if abs(yb_new - yb) < 1e-10, break, end
    if (y2-yb)*(yb_new-yb)>0
        z2 = znew;
    else
        z1 = znew;
    end
end
znew,yb_new-yb
% use znew to run rk4 using dvdtf to get solved BVP
end

function yb_guess = init2bound(z0)
t = 1:0.05:1.4; % time range for simulation (from BC1 to BC2)
y0 = 1; % known variable value at BC1

v0 = [y0 z0]; 
diff_func = @(t,v) dvdtf(t,v);
v = rk4(diff_func,t,v0);
yb_guess = v(end,1); % guessed value for known variable at BC2
end

function dvdt = dvdtf(t,v)
y = v(1); z = v(2);
dvdt(1) = y*t^2 + log(y) + z; % =dy/dt
dvdt(2) = z*t + y; % =dz/dt
end