%% Sliding Mode Control of Inverted Pendulum on Cart

clc; clear all; close all;

%% System Model

%  States: x = [x1, x2]^T = [theta, thetaDot]^T
%
%  | x1 | = |           x2           |
%  | x2 |   | f(x1, x2) + b(x1, x2)u |
%
%  f(x1, x2) = ((M + m) * g * sin(x1) - m * l * x2^2sin(x1)cos(x1)) / ...
%              ((4 / 3) * (M + m) * l - m * g * l * cos^2(x1))
%  b(x1, x2) = cos(x1) / ((4 / 3) * (M + m) * l - m * g * l * cos^2(x1))

%% Controller Design



%% System Specifications

M     = 1; % kg
m     = 0.1; % kg
l     = 0.5; % m
g     = 9.81; % m / s

%% Control System Simulation
x1 = pi; x2 = 0;
x = [pi; 0];
lambda = 5; k = 18;
controlSaturation = 30;

dt = 0.001;
t = dt:dt:10;

x1Des = (pi / 30) .* sin(t);
x2Des = (pi / 30) .* cos(t);
x2DesDot = -(pi / 30) .* sin(t);

for i = 2:length(t)

    e = x1(i - 1) - x1Des(i);
    eDot = x2(i - 1) - x2Des(i);

    s = lambda * e + eDot;
    f = ((M + m) * g * sin(x1(i - 1)) - m * l * x2(i - 1)^2 * sin(x1(i - 1)) * cos(x1(i - 1))) / ((4 / 3) * (M + m) * l - m * l * cos(x1(i - 1))^2);
    b = cos(x1(i - 1)) / ((4 / 3) * (M + m) * l - m * l * cos(x1(i - 1))^2);
    h = lambda * (x2(i - 1) - x2Des(i)) + (f - x2DesDot(i));
    u(:, i) = -(1 / b) * (h + k * sign(s));
    if abs((u(i))) > controlSaturation
        u(i) = sign(u(i)) * controlSaturation;
    end

    xDot1 = x(2, i - 1);
    xDot2 = f + b * u(:, i);
    x1(:, i) = x1(:, i - 1) + xDot1 * dt;
    x2(:, i) = x2(:, i - 1) + xDot2 * dt;
    x(:, i) = [x1(:, i); x2(:, i)];

end


figure()
subplot(2, 1, 1)
plot(t, rad2deg(x1))
hold on
plot(t, rad2deg(x1Des))
legend('Control System', 'Reference')
subplot(2, 1, 2)
plot(t, u)
