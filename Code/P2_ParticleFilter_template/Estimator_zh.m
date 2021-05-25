function [postParticles] = Estimator_zh(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

% Set number of particles:
N_particles = 2000; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    r2 = estConst.d*estConst.d*rand(N_particles,1);
    theta = 2*pi*rand(N_particles,1);
    flag = rand(N_particles,1);
    for i = 1:N_particles
        if (flag(i)<0.5)
            postParticles.x_r(i) = sqrt(r2(i))*cos(theta(i)) + estConst.pA(1); % 1xN_particles matrix
            postParticles.y_r(i) = sqrt(r2(i))*sin(theta(i)) + estConst.pA(2); % 1xN_particles matrix
        else
            postParticles.x_r(i) = sqrt(r2(i))*cos(theta(i)) + estConst.pB(1); % 1xN_particles matrix
            postParticles.y_r(i) = sqrt(r2(i))*sin(theta(i)) + estConst.pB(2); % 1xN_particles matrix
        end
    end
%     if (flag(i)<0.5)
%         postParticles.x_r = sqrt(r2).*cos(theta) + estConst.pA(1); % 1xN_particles matrix
%         postParticles.y_r = sqrt(r2).*sin(theta) + estConst.pA(2); % 1xN_particles matrix
%     else
%         postParticles.x_r= sqrt(r2).*cos(theta) + estConst.pB(1); % 1xN_particles matrix
%         postParticles.y_r = sqrt(r2).*sin(theta) + estConst.pB(2); % 1xN_particles matrix
%     end
    postParticles.phi = 2*estConst.phi_0*rand(N_particles,1) - estConst.phi_0; % 1xN_particles matrix
    postParticles.kappa = 2*estConst.l*rand(N_particles,1) - estConst.l;% 1xN_particles matrix 
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.
% Implement your estimator here!

% Prior Update:
v_f = estConst.sigma_f*rand(N_particles,1) - estConst.sigma_f/2;
v_phi = estConst.sigma_phi*rand(N_particles,1) - estConst.sigma_phi/2;

tempParticles.x_r = prevPostParticles.x_r + (act(1) + v_f) .* cos(prevPostParticles.phi);
tempParticles.y_r = prevPostParticles.y_r + (act(1) + v_f) .* sin(prevPostParticles.phi);
tempParticles.phi = prevPostParticles.phi + act(2) + v_phi;
tempParticles.kappa = prevPostParticles.kappa;


% Posterior Update: 
% weights
% beta = post_weights(N_particles, tempParticles, estConst, sens);
e = estConst.epsilon;
w = zeros(N_particles,1);
beta = ones(N_particles,1) / N_particles;
for i = 1:N_particles
    w(i) = sens - get_distance(tempParticles.x_r(i), tempParticles.y_r(i), tempParticles.phi(i), tempParticles.kappa(i), estConst.contour);
    if (w(i)<(-2.5*e)&&w(i)>=(-3*e))
    beta(i) = (2*w(i)+6*e)/(5*e*e);
    elseif (w(i)<(-2*e)&&w(i)>=(-2.5*e))
    beta(i) = -(4*e+2*w(i))/(5*e*e);
    elseif (w(i)<0&&w(i)>=(-2*e))
    beta(i) = (w(i)+2*e)/(5*e*e);
    elseif (w(i)<(2*e)&&w(i)>=0)
    beta(i) = (-w(i)+2*e)/(5*e*e);
    elseif (w(i)<(2.5*e)&&w(i)>=(2*e))
    beta(i) = (-4*e+2*w(i))/(5*e*e);
    elseif (w(i)<(3*e)&&w(i)>=(2.5*e))
    beta(i) = (-2*w(i)+6*e)/(5*e*e);
    elseif (w(i)<(-3*e)||w(i)>(3*e))
    beta(i) = 1;
    end
end

% if (sum(beta)<0.001) % all particles have zero measurement likelihood
%     for i = 1:N_particles
%         x = 0.04*rand(1,1)-0.02;
%         y = 0.04*rand(1,1)-0.02;
%         tempParticles.x_r(i) = tempParticles.x_r(i) + x;
%         tempParticles.y_r(i) = tempParticles.y_r(i) + y;
%         tempParticles.phi(i) = 0.04*pi*rand(1,1)-0.02*pi;
%         tempParticles.kappa(i) = 0.1*estConst.l*rand(1,1)-0.05*estConst.l;
%         if (tempParticles.x_r(i)>3)
%             tempParticles.x_r(i)=3;
%         end
%         if (tempParticles.x_r(i)<-estConst.l)
%             tempParticles.x_r(i) = -estConst.l;
%         end
%         if (tempParticles.y_r(i)>3)
%             tempParticles.y_r(i) = 3;
%         end
%         if (tempParticles.y_r(i)<0)
%             tempParticles.y_r(i)=0;
%         end
%         if(tempParticles.kappa(i)<-estConst.l)
%             tempParticles.kappa(i)=-estConst.l;
%         end
%         if(tempParticles.kappa(i)>estConst.l)
%             tempParticles.kappa(i)=estConst.l;
%         end
%         beta(i) = 1/N_particles;
%     end
% end

beta = beta / sum(beta);%normalize


%resample
for i = 1:N_particles
    r = rand;
    beta_sum = 0;
    for m = 1:N_particles
        beta_sum = beta_sum + beta(m);
        if (beta_sum>=r)
            postParticles.x_r(i) = tempParticles.x_r(m);
            postParticles.y_r(i) = tempParticles.y_r(m);
            postParticles.phi(i) = tempParticles.phi(m);
            postParticles.kappa(i) = tempParticles.kappa(m);
            break;
        end
    end
end
%roughening
k = 0.0001;
d = 4;
delta1 = k*N_particles^(-1/d)*(max(postParticles.x_r)-min(postParticles.x_r));
delta2 = k*N_particles^(-1/d)*(max(postParticles.y_r)-min(postParticles.y_r));
delta3 = k*N_particles^(-1/d)*(max(postParticles.phi)-min(postParticles.phi));
delta4 = k*N_particles^(-1/d)*(max(postParticles.kappa)-min(postParticles.kappa));

postParticles.x_r = postParticles.x_r + delta1*randn(N_particles,1)';
postParticles.y_r = postParticles.y_r + delta2*randn(N_particles,1)';
postParticles.phi = postParticles.phi + delta3*randn(N_particles,1)';
postParticles.kappa = postParticles.kappa + delta4*randn(N_particles,1)';
end % end estimator


function [beta] = post_weights(N_particles, tempParticles, estConst, sens)
e = estConst.epsilon;
w = zeros(N_particles,1);
beta = ones(N_particles,1)/N_particles;
for i = 1:N_particles
    w(i) = sens - get_distance(tempParticles.x_r(i), tempParticles.y_r(i), tempParticles.phi(i), tempParticles.kappa(i), estConst.contour);
    if (w(i)<(-2.5*e)&&w(i)>=(-3*e))
    beta = (2*w(i)+6*e)/(5*e*e);
    elseif (w(i)<(-2*e)&&w(i)>=(-2.5*e))
    beta = -(4*e+2*w(i))/(5*e*e);
    elseif (w(i)<0&&w(i)>=(-2*e))
    beta = (w(i)+2*e)/(5*e*e);
    elseif (w(i)<(2*e)&&w(i)>=0)
    beta = (-w(i)+2*e)/(5*e*e);
    elseif (w(i)<(2.5*e)&&w(i)>=(2*e))
    beta = (-4*e+2*w(i))/(5*e*e);
    elseif (w(i)<(3*e)&&w(i)>=(2.5*e))
    beta = (-2*w(i)+6*e)/(5*e*e);
    elseif (w(i)<(-3*e)||w(i)>(3*e))
    beta = 0;
    end
end
end

%% get_distance
function [dis] = get_distance(xr, yr, phi, kappa, contour)
    xc = 100;
    yc = 100;
    contour(8,1) = kappa;
    contour(9,1) = kappa;
    contour = [contour; contour(1,:)]; % add first point in contour to make convenience for the following loop
    for i = 1:10
        [x, y] = inter_line([xr, xr+5*cos(phi)], [yr, yr+5*sin(phi)], [contour(i,1),contour(i+1,1)],[contour(i,2),contour(i+1,2)]);
        if x~=0 || y~=0
            if (x-xr)^2+(y-yr)^2 < (x-xc)^2+(y-yc)^2
                xc = x;
                yc = y;
            end
        end
    end
    dis = sqrt((xc-xr)^2+(yc-yr)^2);
end

%% interline
function [x,y] = inter_line(x1, y1, x2, y2)
    if((x1(1)~=x1(2))&&(x2(1)~=x2(2)))
        k1 = (y1(1)-y1(2))/(x1(1)-x1(2));
        k2 = (y2(1)-y2(2))/(x2(1)-x2(2));
        if(abs(k1-k2)< 0.01)
            x=0;
            y=0;
        else
            x = (y2(2)-y1(2)+k1*x1(2)-k2*x2(2))/(k1-k2);
            y = k1*x-k1*x1(2)+y1(2);
        end
    elseif((x1(1)==x1(2))&&(x2(1)~=x2(2)))
        x = x1(1);
        y = (y2(1)-y2(2))*(x-x2(2))/(x2(1)-x2(2))+y2(2);
    elseif((x1(1)~=x1(2))&&(x2(1)==x2(2)))
        x = x2(1);
        y = (y1(1)-y1(2))*(x-x1(2))/(x1(1)-x1(2))+y1(2);
    else
        x = 0;
        y = 0;
    end
    if(((x-x1(1))*(x-x1(2))>0 && x1(1)~=x1(2))||((x-x2(1))*(x-x2(2))>0&&x2(1)~=x2(2))||((y-y1(1))*(y-y1(2))>0&&y1(1)~=y1(2))||((y-y2(1))*(y-y2(2))>0&&y2(1)~=y2(2)))
        x=0;
        y=0;
    end
end


%% intersection
function [x_intersec, y_intersec] = intersection(x1, y1, x2, y2)
L1_x1 = x1(1);
L1_x2 = x1(2);
L1_y1 = y1(1);
L1_y2 = y1(2);

L2_x1 = x2(1);
L2_x2 = x2(2);
L2_y1 = y2(1);
L2_y2 = y2(2);

A = [L1_x2 - L1_x1, -(L2_x2 - L2_x1);...
    L1_y2 - L1_y1, -(L2_y2 - L2_y1)];
b = [L2_x1 - L1_x1; L2_y1 - L1_y1];
t = A \ b;
x_intersec = L1_x1 + t(1) * (L1_x2 - L1_x1);
y_intersec = L1_y1 + t(1) * (L1_y2 - L1_y1);
end

% function [dis] = get_distance(xr, yr, phi, kappa, contour)
%     x = [xr, xr+5*cos(phi)];
%     y = [yr, yr+5*sin(phi)];
%     polyx = contour(:,1);
%     polyy = contour(:,2);
%     polyx(8) =  kappa;
%     polyx(9) =  kappa;
%     polyx = [polyx;polyx(1)];
%     polyy = [polyy;polyy(1)];
%     [xi, yi] = polyxpoly(x,y,polyx,polyy);
%     dis = 100;
%     for i = 1:size(xi,1)
%         if(sqrt((xi(i)-xr)*(xi(i)-xr)+(yi(i)-yr)*(yi(i)-yr))<dis)
%             dis = sqrt((xi(i)-xr)*(xi(i)-xr)+(yi(i)-yr)*(yi(i)-yr));
%         end
%     end
% end

