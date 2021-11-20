function IVP_test
t=1:0.1:1.4;
f0 = [1];
df = @(t,f) fun(t,f);
f_euler = euler_IVP(df,t,f0);
f_rk2 = rk2(df,t,f0);
t2=1:0.2:1.4;
f_rk4 = rk4(df,t2,f0);

stiff_demo
end

% Sample function 'fun'
function dfdt = fun(t,f)
dfdt(1) = f(1)*t^2;
end

function stiff_demo
h = 0.1; start = 0; endpoint = 1;
t = start:h:endpoint;
y0 = 2; z0 = 1;
[t,y,z] = IVP_implicit_example(t,y0,z0);
subplot(1,2,1)
    plot(t,[y z])
    xlabel('t')
    ylabel('value of y and z')
    legend('y','z')
    title('implicit method')

[t,y,z] = IVP_explicit_example(t,y0,z0);
subplot(1,2,2)
    plot(t,[y z])
    xlabel('t')
    ylabel('value of y and z')
    legend('y','z')
    title('explicit method')
    
end

function [t,y,z] = IVP_explicit_example(t,y0,z0)
h = t(2) - t(1);
a = -500.5;
b = 499.5;
c = 499.5;
d = -500.5;
k1 = length(t);
dfdt = @(t,y,z) a*y+b*z;
dzdt = @(t,y,z) c*y+d*z;
y = y0; z = z0;
for i=1:k1-1
    t1 = t(i); y1 = y(i); z1 = z(i);
    y2 = y1 + h*dfdt(t1,y1,z1);
    z2 = z1 + h*dzdt(t1,y1,z1);
    y(i+1,1) = y2;
    z(i+1,1) = z2;
end
end

