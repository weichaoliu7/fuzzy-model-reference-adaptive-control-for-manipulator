function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 2*3*5;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 10;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = zeros(1,2*3*5);
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
q1_desired = u(1); % desired angular position of joint 1
dq1_desired = u(2); % desired angular velocity of joint 1
q2_desired = u(3); % desired angular position of joint 2
dq2_desired = u(4); % desired angular velocity of joint 2

q1 = u(5); % actual angular position of joint 1
dq1 = u(6); % actual angular velocity of joint 1
q2 = u(7); % actual angular position of joint 2
dq2 = u(8); % actual angular velocity of joint 2

xi1 = [q1;dq1]; % state variable of the first reference model
xi2 = [q2;dq2]; % state variable of the second reference model
xi = [xi1;xi2]; % state variable

r1 = q1_desired; % bounded reference input of the first reference model
r2 = q2_desired; % bounded reference input of the second reference model

z1 = [xi;r1]; % defined in Eq. 5
z2 = [xi;r2];

error1 = [r1 - xi1(1);dq1_desired-xi1(2)]; % tracking error of the first reference model
error2 = [r2 - xi2(1);dq2_desired-xi2(2)]; % tracking error of the second reference model

% membership function of input vectors in fuzzy sets
membership = struct('q1', q1, 'q2', q2);

% membership function of q1
% first membership function of q1, V11, implies that q1 is small
if xi(1)<=0 && xi(1)>=-pi/2
    membership.q1(1) = - xi(1)/(pi/2);
else
    membership.q1(1) = 0;
end

% second membership function of q1, V12, implies that q1 is medium
if xi(1)<=0 && xi(1)>=-pi/2
    membership.q1(2) = (xi(1)+pi/2)/(pi/2);
elseif xi(1)<=pi && xi(1)>=0
    membership.q1(2) = (pi/2-xi(1))/(pi/2);    
else
    membership.q1(2) = 0;
end

% third membership function of q1, V13, implies that q1 is large
if xi(1)<=pi && xi(1)>=0
    membership.q1(3) = xi(1)/(pi/2);
else
    membership.q1(3) = 0;
end

% membership function of q2
% first membership function of q2, V21, implies that q2 is small
if xi(3)<=0 && xi(3)>=-pi/2
    membership.q2(1) = - xi(3)/(pi/2);
else
    membership.q2(1) = 0;
end

% second membership function of q2, V22, implies that q2 is medium
if xi(3)<=0 && xi(3)>=-pi/2
    membership.q2(2) = (xi(3)+pi/2)/(pi/2);
elseif xi(3)<=pi && xi(3)>=0
    membership.q2(2) = (pi/2-xi(3))/(pi/2);
else
    membership.q2(2) = 0;
end

% third membership function of q2, V23, implies that q2 is large
if xi(3)<=pi && xi(3)>=0
    membership.q2(3) = xi(3)/(pi/2);
else
    membership.q2(3) = 0;
end

sum = membership.q1(1)+membership.q1(2)+membership.q1(3);
xi_1 = [membership.q1(1) membership.q1(2) membership.q1(3)]./sum; % normalized firing strength vector, defined in Eq. 6
sum = membership.q2(1)+membership.q2(2)+membership.q2(3);
xi_2 = [membership.q2(1) membership.q2(2) membership.q2(3)]./sum;

gamma12 = 20; % fuzzy controller parameter update gain, 7
gamma22 = 1000; % 300
bc1 = [0;1];
bc2 = [0;1];

Am1 = [0 1;-16 -8]; % state matrix of the first reference model
Am2 = [0 1;-16 -8]; % state matrix of the second reference model
Q1 = [15 0;0 5];
Q2 = [15 0;0 5];

% solution of reference model Lyapunov equation
P1 = lyap(Am1', Q1);
P2 = lyap(Am2', Q2);

theta1_derivative = gamma12 * bc1' * P1 * error1 * xi_1' * z1'; % Eq. 9
theta2_derivative = gamma22 * bc2' * P2 * error2 * xi_2' * z2';

% handle non-numeric floating point value
for i=1:3
    for j=1:5
        if isnan(theta1_derivative(i,j))
            theta1_derivative(i,j)=0;
        end
        if isnan(theta2_derivative(i,j))
            theta2_derivative(i,j)=0;
        end
    end
end

for i=1:3
    for j=1:5
        sys(5*(i-1)+j) = theta1_derivative(i,j);
        sys(15+5*(i-1)+j) = theta2_derivative(i,j);
    end
end

function sys=mdlOutputs(t,x,u)
q1_desired = u(1); % desired angular position of joint 1
dq1_desired = u(2); % desired angular velocity of joint 1
q2_desired = u(3); % desired angular position of joint 2
dq2_desired = u(4); % desired angular velocity of joint 2

q1 = u(5); % actual angular position of joint 1
dq1 = u(6); % actual angular velocity of joint 1
q2 = u(7); % actual angular position of joint 2
dq2 = u(8); % actual angular velocity of joint 2

torque1=u(9); % control input of the first reference model
torque2=u(10); % control input of the second reference model

xi1 = [q1;dq1]; % state variable of the first reference model
xi2 = [q2;dq2]; % state variable of the second reference model
xi = [xi1;xi2]; % state variable

r1 = q1_desired;% bounded reference input of the first reference model
r2 = q2_desired;% bounded reference input of the second reference model

z1 = [xi;r1]; % defined in Eq. 5
z2 = [xi;r2];

error1 = [r1 - xi1(1);dq1_desired-dq1]; % tracking error of the first reference model
error2 = [r2 - xi2(1);dq2_desired-dq2]; % tracking error of the second reference model

% membership function of input vectors in fuzzy sets
membership = struct('q1', q1, 'q2', q2);

% membership function of q1
% first membership function of q1, V11, implies that q1 is small
if xi(1)<=0 && xi(1)>=-pi/2
    membership.q1(1) = - xi(1)/(pi/2);
else
    membership.q1(1) = 0;
end

% second membership function of q1, V12, implies that q1 is medium
if xi(1)<=0 && xi(1)>=-pi/2
    membership.q1(2) = (xi(1)+pi/2)/(pi/2);
elseif xi(1)<=pi && xi(1)>=0
    membership.q1(2) = (pi/2-xi(1))/(pi/2);    
else
    membership.q1(2) = 0;
end

% third membership function of q1, V13, implies that q1 is large
if xi(1)<=pi && xi(1)>=0
    membership.q1(3) = xi(1)/(pi/2);
else
    membership.q1(3) = 0;
end

% membership function of q2
% first membership function of q2, V21, implies that q2 is small
if xi(3)<=0 && xi(3)>=-pi/2
    membership.q2(1) = - xi(3)/(pi/2);
else
    membership.q2(1) = 0;
end

% second membership function of q2, V22, implies that q2 is medium
if xi(3)<=0 && xi(3)>=-pi/2
    membership.q2(2) = (xi(3)+pi/2)/(pi/2);
elseif xi(3)<=pi && xi(3)>=0
    membership.q2(2) = (pi/2-xi(3))/(pi/2);
else
    membership.q2(2) = 0;
end

% third membership function of q2, V23, implies that q2 is large
if xi(3)<=pi && xi(3)>=0
    membership.q2(3) = xi(3)/(pi/2);
else
    membership.q2(3) = 0;
end
sum = membership.q1(1)+membership.q1(2)+membership.q1(3);
xi_1 = [membership.q1(1) membership.q1(2) membership.q1(3)]./sum; % normalized firing strength vector
sum = membership.q2(1)+membership.q2(2)+membership.q2(3);
xi_2 = [membership.q2(1) membership.q2(2) membership.q2(3)]./sum;

gamma11 = 100; % update gain, 7
gamma21 = 10; % 300

bc1 = [0;1];
bc2 = [0;1];

Am1 = [0 1;-16 -8]; % state matrix of the first reference model
Am2 = [0 1;-16 -8]; % state matrix of the second reference model
Q1 = [15 0;0 5];
Q2 = [15 0;0 5];

% solution of reference model Lyapunov equation
P1 = lyap(Am1', Q1);
P2 = lyap(Am2', Q2);

% fuzzy controller parameter update PI law, Eq. 37
phi1 = gamma11 * bc1' * P1 * error1 * xi_1' * z1';% proportional term, Eq. 8
phi2 = gamma21 * bc2' * P2 * error2 * xi_2' * z2';

theta1 = zeros(3,5);
theta2 = zeros(3,5);

for i=1:3
    for j=1:5
        theta1(i,j) = x(5*(i-1)+j); % integral term of parameter update PI law
        theta2(i,j) = x(15+5*(i-1)+j);
    end
end

THETA1 = phi1 + theta1;
THETA2 = phi2 + theta2;

% handle of fuzzy controller parameter out-of-range
for i=1:3
    for j=1:5
        if THETA1(i,j) >= 60
            THETA1(i,j) = 60;
        elseif THETA1(i,j) <= -60
            THETA1(i,j) = -60;
        end
        if THETA2(i,j) >= 60
            THETA2(i,j) = 60;
        elseif THETA2(i,j) <= -60
            THETA2(i,j) = -60;
        end
    end
end

% handle non-numeric floating point value of fuzzy controller parameter
for i=1:3
    for j=1:5
        if isnan(THETA1(i,j))
            THETA1(i,j) = 0;
        end
        if isnan(THETA1(i,j))
            THETA1(i,j) = 0;
        end
    end
end

% fuzzy logic system obtained through the combination of single-point fuzzification, product inference, and center-weighted defuzzification
torque1_fuzzy = xi_1 * THETA1 * z1; % output of the first fuzzy controller, Eq. 50
torque2_fuzzy = xi_2 * THETA2 * z2; % output of the second fuzzy controller, Eq. 51

% handle non-numeric floating point value of output of fuzzy controller
if isnan(torque1_fuzzy)
    torque1_fuzzy = u(9);
end

if isnan(torque2_fuzzy)
    torque1_fuzzy = u(10);
end

b11_lower = 1.7; % lower bound of bii,smooth unknown functions
b22_lower = 2;
b12_upper = 1.3; % lower bound of bij,i != j,smooth unknown functions
b21_upper = 1.3;

omega1_upper = 0.9; % upper bound of minimum approximation error
omega2_upper = 0.9;

eta1_upper = 0.78; % upper bound of external disturbance
eta2_upper = 3.13;

d1_upper = omega1_upper + eta1_upper / b11_lower; % upper bound of uncertain term
d2_upper = omega2_upper + eta2_upper / b22_lower;

delta1 = 0.1;
delta2 = 0.1;
beta1 = 0.1;
beta2 = 0.1;

Beta11 = 0.8 * abs(xi(4)); % known function
Beta22 = 6.5 * abs(xi(4));

% additional control term used to overcome uncertainties
temp11 = (1/b11_lower * (b12_upper * abs(torque2) + d1_upper)) * (1+(delta1/beta1));
temp12 = ((error1'*P1*bc1)/(abs(error1'*P1*bc1)+delta1));
temp13 = Beta11/(2*b11_lower^2) * error1' * P1 * error1;
torque1_si = temp11 * temp12 + temp13;

temp21 = (1/b22_lower * (b21_upper * abs(torque1) + d2_upper)) * (1+(delta2/beta2));
temp22 = ((error2'*P2*bc2)/(abs(error2'*P2*bc2)+delta2));
temp23 = Beta22/(2*b22_lower^2) * error2' * P2 * error2;
torque2_si = temp21 * temp22 + temp23;

% control input torque, Eq. 11
torque1 = torque1_fuzzy + torque1_si;
torque2 = torque2_fuzzy + torque2_si;

% handle of control input torque out-of-range
if torque1 >= 45
    torque1 = u(9);
elseif torque1 <= -45
    torque1 = u(9);
end
if torque2 >= 45
    torque2 = u(10);
elseif torque2 <= -45
    torque2 = u(10);
end

% handle non-numeric floating point value of control input torque
if isnan(torque1)
    torque1 = u(9);
end

if isnan(torque2)
    torque2 = u(10);
end

sys(1)=torque1;
sys(2)=torque2;