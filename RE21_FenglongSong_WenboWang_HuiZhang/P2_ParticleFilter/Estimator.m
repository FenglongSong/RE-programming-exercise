function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
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
N_particles = 5000; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    r2 = estConst.d*estConst.d*rand(1,N_particles);
    theta = 2*pi*rand(1,N_particles);
    flag = rand(1,N_particles);
    for i = 1:N_particles
        if (flag(i)<0.5)
            postParticles.x_r(i) = sqrt(r2(i))*cos(theta(i)) + estConst.pA(1); % 1xN_particles matrix
            postParticles.y_r(i) = sqrt(r2(i))*sin(theta(i)) + estConst.pA(2); % 1xN_particles matrix
        else
            postParticles.x_r(i) = sqrt(r2(i))*cos(theta(i)) + estConst.pB(1); % 1xN_particles matrix
            postParticles.y_r(i) = sqrt(r2(i))*sin(theta(i)) + estConst.pB(2); % 1xN_particles matrix
        end
    end
    postParticles.phi = 2*estConst.phi_0*rand(1,N_particles)-estConst.phi_0*ones(1,N_particles); % 1xN_particles matrix
    postParticles.kappa = 2*estConst.l*rand(1,N_particles)-estConst.l*ones(1,N_particles);% 1xN_particles matrix 
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.
% Implement your estimator here!
% Prior Update:
low_variance_flag = 0; % 0:normal reample   1:low variance resample
v_f = estConst.sigma_f*rand(1,N_particles)-estConst.sigma_f/2*ones(1,N_particles);
v_phi = estConst.sigma_phi*rand(1,N_particles)-estConst.sigma_phi/2*ones(1,N_particles);
for i = 1:N_particles
    tempParticles.x_r(i) = prevPostParticles.x_r(i) + (act(1) + v_f(i)) * cos(prevPostParticles.phi(i));
    tempParticles.y_r(i) = prevPostParticles.y_r(i) + (act(1) + v_f(i)) * sin(prevPostParticles.phi(i));
    tempParticles.phi(i) = prevPostParticles.phi(i) + act(2) + v_phi(i);
    tempParticles.kappa(i) = prevPostParticles.kappa(i);
end
% Posterior Update: 
% weights
[beta] = post_weights(N_particles, tempParticles, estConst, sens);

if (sum(beta)==0) % all particles have zero measurement likelihood
    disp("lost!");
    tempParticles.x_r = tempParticles.x_r+0.2*randn(1,N_particles)-0.1*ones(1,N_particles);
    tempParticles.y_r = tempParticles.y_r+0.2*randn(1,N_particles)-0.1*ones(1,N_particles);
    tempParticles.phi = tempParticles.phi+0.2*pi*rand(1,N_particles)-0.1*pi*ones(1,N_particles);
    tempParticles.kappa = tempParticles.kappa+0.6*estConst.l*rand(1,N_particles)-0.3*estConst.l*ones(1,N_particles);
    [tempParticles, beta] = put_inside(tempParticles, N_particles, estConst, beta);
    
    [beta] = post_weights(N_particles, tempParticles, estConst, sens);
    if (sum(beta)==0)
        disp("still lost!");
        low_variance_flag = 1;
        tempParticles.x_r = tempParticles.x_r+0.4*randn(1,N_particles)-0.2*ones(1,N_particles);
        tempParticles.y_r = tempParticles.y_r+0.4*randn(1,N_particles)-0.2*ones(1,N_particles);
        tempParticles.phi = tempParticles.phi+0.6*pi*rand(1,N_particles)-0.3*pi*ones(1,N_particles);
        tempParticles.kappa = tempParticles.kappa+0.6*estConst.l*rand(1,N_particles)-0.3*estConst.l*ones(1,N_particles);
        
        [tempParticles, beta] = put_inside(tempParticles, N_particles, estConst, beta);
    end
    [beta] = post_weights(N_particles, tempParticles, estConst, sens);
    if (sum(beta)==0)
        beta = ones(1,N_particles);
    end    
end
beta = beta/sum(beta);%normalize

if low_variance_flag
%low variance resample (it's faster)
    r = rand * 1/N_particles;                       % Select random number between 0-1
    w = beta(1);                  % Initial weight
    i = 1;
    j = 1;

    for m = 1:N_particles
        U = r + (m - 1)/N_particles; % Index of original sample + size^-1
        while U > w                 % I'm not sure what this loop is doing
            i = i + 1;
            w = w + beta(i);
        end
        postParticles.x_r(j) = tempParticles.x_r(i);   % Add selected sample to resampled array
        postParticles.y_r(j) = tempParticles.y_r(i);
        postParticles.phi(j) = tempParticles.phi(i);
        postParticles.kappa(j) = tempParticles.kappa(i);
        j = j + 1;
    end
else
%normal resample
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
end
%roughening
k = 0.1;
d = 4;
delta1 = k*N_particles^(-1/d)*(max(postParticles.x_r)-min(postParticles.x_r));
delta2 = k*N_particles^(-1/d)*(max(postParticles.y_r)-min(postParticles.y_r));
delta3 = k*N_particles^(-1/d)*(max(postParticles.phi)-min(postParticles.phi));
delta4 = k*N_particles^(-1/d)*(max(postParticles.kappa)-min(postParticles.kappa));

postParticles.x_r = postParticles.x_r + delta1*randn(N_particles,1)';
postParticles.y_r = postParticles.y_r + delta2*randn(N_particles,1)';
postParticles.phi = postParticles.phi + delta3*randn(N_particles,1)';
postParticles.kappa = postParticles.kappa + delta4*randn(N_particles,1)';
if(postParticles.kappa(i)<-estConst.l)
    postParticles.kappa(i)=-estConst.l;
end
if(postParticles.kappa(i)>estConst.l)
    postParticles.kappa(i)=estConst.l;
end
end % end estimator

function [tempParticles, beta] = put_inside(tempParticles, N_particles, estConst, beta)
for i = 1:N_particles
    if (tempParticles.x_r(i)>3)
    tempParticles.x_r(i)=3;
    end
    if (tempParticles.x_r(i)<-estConst.l)
        tempParticles.x_r(i) = -estConst.l;
    end
    if (tempParticles.y_r(i)>3)
        tempParticles.y_r(i) = 3;
    end
    if (tempParticles.y_r(i)<0)
        tempParticles.y_r(i)=0;
    end
    if(tempParticles.kappa(i)<-estConst.l)
        tempParticles.kappa(i)=-estConst.l;
    end
    if(tempParticles.kappa(i)>estConst.l)
        tempParticles.kappa(i)=estConst.l;
    end
    if ((tempParticles.x_r(i)+tempParticles.y_r(i))>5)
        beta(i)=0;
    end
end
end

function [dis] = get_distance(xr, yr, phi, kappa, contour)
    xc = 100;
    yc = 100;
    contour(8,1) = kappa;
    contour(9,1) = kappa;
    contour = [contour;contour(1,:)];
    for i=1:10
        [x,y] = inter_line([xr, xr+5*cos(phi)], [yr, yr+5*sin(phi)], [contour(i,1),contour(i+1,1)],[contour(i,2),contour(i+1,2)]);
        if(x~=0||y~=0)
            if(((x-xr)^2+(y-yr)^2)<((x-xc)^2+(y-yc)^2))
                xc = x;
                yc = y;
            end
        end
    end
    dis = sqrt((xc-xr)^2+(yc-yr)^2);
end

function [x,y] = inter_line(x1,y1,x2,y2)
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



function [beta] = post_weights(N_particles, tempParticles, estConst, sens)
e = estConst.epsilon;
w = zeros(N_particles,1);
beta = ones(N_particles,1)/N_particles;
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
    beta(i) = 0;
    end
end
end