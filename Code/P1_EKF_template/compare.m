% compare
px = 1; py = 2;
sx = 1; sy = 2;
phi = 5; rho = 6; b = 10;
u_t = 45; u_r = 65;

estConst = EstimatorConst();
Cdh = estConst.dragCoefficientHydr;
Cda = estConst.dragCoefficientAir;
Cr = estConst.rudderCoefficient;
Cw = estConst.windVel;

A = zeros(7);
squareroot = sqrt( (sx-Cw*cos(rho))^2 + (sy-Cw*sin(rho))^2 );
A(1,3) = 1;
A(2,4) = 1;
A(3,3) = cos(phi) * (-2 * Cdh * sx) - Cda * squareroot - Cda * (sx-Cw*cos(rho))^2 / squareroot;
A(3,4) = cos(phi) * (-2 * Cdh * sy) - Cda * (sx-Cw*cos(rho))*(sy-Cw*sin(rho)) / squareroot;
A(3,5) = -sin(phi) * (tanh(u_t) - Cdh*(sx^2+sy^2));
A(3,6) = -Cda * Cw * sin(rho) * squareroot + Cda * (Cw*cos(rho) - sx) * Cw*(sx*sin(rho) - sy*cos(rho)) / squareroot;
A(4,3) = sin(phi) * (-2 * Cdh * sx) - Cda * (sx-Cw*cos(rho))*(sy-Cw*sin(rho)) / squareroot;
A(4,4) = sin(phi) * (-2 * Cdh * sy) - Cda * squareroot - Cda * (sy-Cw*sin(rho))^2 / squareroot;
A(4,5) = cos(phi) * (tanh(u_t) - Cdh*(sx^2+sy^2));
A(4,6) = Cda * Cw * cos(rho) * squareroot + Cda * (Cw*sin(rho) - sy) * Cw*(sx*sin(rho) - sy*cos(rho)) / squareroot;

B = A;

A = zeros(7, 7);
A(1, 3) = 1;
A(2, 4) = 1;
A(3, 3) = -2*cos(phi)*Cdh*sx - Cda*(squareroot + (sx-Cw*cos(rho))*(sx-Cw*cos(rho))/squareroot);
A(4, 3) = -2*sin(phi)*Cdh*sx - Cda*(squareroot + (sx-Cw*sin(rho))*(sx-Cw*cos(rho))/squareroot);
A(3, 4) = -2*cos(phi)*Cdh*sy - Cda*(sx-Cw*cos(rho))*(sy-Cw*sin(rho))/squareroot;
A(4, 4) = -2*sin(phi)*Cdh*sy - Cda*(sx-Cw*cos(rho))*(sy-Cw*sin(rho))/squareroot;
A(3, 5) = -sin(phi)*tanh(u_t) + sin(phi)*Cdh*(sx*sx+sy*sy);
A(4, 5) = cos(phi)*tanh(u_t) - cos(phi)*Cdh*(sx*sx+sy*sy);
A(3, 6) = -Cda*Cw*sin(phi)*squareroot - Cda*(sx-Cw*cos(rho))*((sx-Cw*cos(rho))*Cw*sin(rho)-(sy-Cw*sin(rho))*Cw.*cos(rho))/squareroot;
A(4, 6) = Cda*Cw*cos(phi)*squareroot - Cda*(sx-Cw*sin(rho))*((sx-Cw*cos(rho))*Cw*sin(rho)-(sy-Cw*sin(rho))*Cw*cos(rho))/squareroot;

A-B