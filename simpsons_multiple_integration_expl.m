function simpsons_multiple_integration_expl
 
f = @(x,y) 2 - (x.^2 + y.^2);
 
xdata = -1:0.1:1;
ydata = -1:0.1:1;
 
k = length(xdata);
 
for i=1:k
    xi = xdata(i);
    zdatai = f(xi, ydata);
    sectionarea(i) = simpsons(ydata, zdatai);
end
 
volume = simpsons(xdata, sectionarea)
 
end