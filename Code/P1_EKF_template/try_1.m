close all;
x = linspace(-1, 2, 200);
y = linspace(-1, 1, 200);
[X,Y]=meshgrid(x,y);
Z=X.^2 + Y.^2 + 10.* (X-1)^2;%
surf(X,Y,Z);
title('{$f(x,y)=x^{2}+y^{2}$}','interpreter','latex');
xlabel('x');
ylabel('y');