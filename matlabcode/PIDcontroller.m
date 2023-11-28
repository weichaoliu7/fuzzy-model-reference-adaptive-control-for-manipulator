function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
switch flag
case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
% case 1
%     sys=mdlDerivatives(t,x,u);
case 3
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 8;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [];

function sys=mdlOutputs(t,x,u)
q1_desired = u(1); % desired angular position of joint 1
dq1_desired = u(2); % desired angular velocity of joint 1
q2_desired = u(3); % desired angular position of joint 2
dq2_desired = u(4); % desired angular velocity of joint 2

q1 = u(5); % actualangular position of joint 1
dq1 = u(6); % actual angular velocity of joint 1
q2 = u(7); % actual angular position of joint 2
dq2 = u(8); % actual angular velocity of joint 2

error1_d = q1 - q1_desired; % position tracking error of joint 1
error2_d = q2 - q2_desired; % position tracking error of joint 2
error1_derivative_d = dq1 - dq1_desired; % velocity tracking error of joint 1
error2_derivative_d = dq2 - dq2_desired; % velocity tracking error of joint 2
error_d = [error1_d;error2_d];
error_derivative_d = [error1_derivative_d;error2_derivative_d];
Kd = [100,0;0,100];
Kp = [2000,0;0,2000];

% control input torque
torque_r =  -(Kp*error_d + Kd*error_derivative_d); % PD controller

sys(1)=torque_r(1);
sys(2)=torque_r(2);