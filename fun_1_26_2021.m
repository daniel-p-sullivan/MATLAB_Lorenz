addpath('/usr/local/MATLAB/R2019b/toolbox/matlab/graphics/');

%% system parameters

rho = [14 13 15 28];
sigma = 10;
beta = 8/3;

%% time simulation parameters
dt = 0.01;
end_time = 100;
t_span = 0:dt:end_time;

%% initial conditions

c0_1 = [1; 1; 1];
c0_2 = [2; 2; 2];


for i = 1:length(rho)
    %% ode solver

    [t, c1] = ode45(@(t, c) lorenz(t, c, rho(i), sigma, beta), t_span, c0_1);
    [t, c2] = ode45(@(t, c) lorenz(t, c, rho(i), sigma, beta), t_span, c0_2);

    %% find fixed points
    
    %fixed_point = fsolve(@(d) lorenz2(d, rho(i), sigma, beta), c0);
    
    %% plot
    subplot(2, 2, i)
    sgtitle('Solutions for the Lorenz System')
    hold on
    plot3(c1(:, 1), c1(:, 2), c1(:, 3), 'blue')
    plot3(c2(:, 1), c2(:, 2), c2(:, 3), 'red')
    plot3(c0_1(1), c0_1(2), c0_1(3), 'black*')
    plot3(c0_2(1), c0_2(2), c0_2(3), 'black*')
    title(['\rho = ' num2str(rho(i)) ', \sigma = ' num2str(sigma) ', \beta = ' sprintf('%.3f', beta)])

end




%% strange attractor equation

function dy_dt = lorenz(t, c, rho, sigma, beta)
    dy_dt = [sigma * (c(2) - c(1));
             c(1)*(rho - c(3))-c(2);
             c(1)*c(2) - beta*c(3)];
end

function dy_dt = lorenz2(c, rho, sigma, beta)
    dy_dt = [sigma * (c(2) - c(1));
             c(1)*(rho - c(3))-c(2);
             c(1)*c(2) - beta*c(3)];

end

