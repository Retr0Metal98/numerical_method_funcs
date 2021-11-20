function [t,y,z] = IVP_implicit_example(t,y0,z0)
%IVP_IMPLICIT_EXAMPLE example of implicit method for stiff IVPs 
k = length(t);
h = t(2) - t(1);
% IVPs: extend with more equations if needed
a = -500.5; b = 499.5; % y' = a*y + b*z
c = 499.5; d = -500.5; % z' = c*y + d*z

% implicit method (see Lecture 11, slide 43-45)
mat = [1-a*h -b*h;
       -c*h 1-d*h;];
y = y0; z = z0; % initialize solution vectors with init values.
for i=1:k-1
    y1 = y(i); z1 = z(i);
    vec = [y1 z1]';
    R_aug = rref_manual([mat vec]);
    sol = R_aug(:,end);
    y2 = sol(1); z2 = sol(2);
    y(i+1,1) = y2;
    z(i+1,1) = z2;
end
end

