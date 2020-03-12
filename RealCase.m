% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

clc; clear; close 
addpath('./utils');

%% load time-series data
% from website: http://data.stats.gov.cn/easyquery.htm?cn=C01
xtrain = flip([ 64.50 57.46 49.94 44.60 44.80...
                33.12 31.16 28.74 25.70 22.70 21.96 17.20]'); 
xtest = flip([nan nan 86.40 81.20 70.85]'); 

ttrain = (1:length(xtrain))';               % time index 
ttest = length(xtrain)+(1:length(xtest))';

%% integral matching-based ode model
poly_order = 1;                 % d/dt x = a*x+b3*t^3+b2*t^2+b1*t^1+c

% step 1: data preprocessing 
[nobs,nsta] = size(xtrain);
tgap = diff(ttrain); 
Theta = tgap.*ones(1,nsta).*cumsum(xtrain(1:nobs-1,:)+xtrain(2:nobs,:))/2;
Const = ones(nobs-1,1);

% step 2: parameter estimation d/dt x = a*x+b3*t^3+b2*t^2+b1*t^1+c
switch poly_order
    case -1     % d/dt x = a*x
        pars = regress(xtrain(2:nobs), ...
                    [Theta, ...                                     % state term 
                    Const]);                                        % initial condition 
    case 0      % d/dt x = a*x+c
        pars = regress(xtrain(2:nobs), ...
                    [Theta, ...                                     % state term 
                    Const, ...                                      % initial condition 
                    ttrain(2:nobs)-ttrain(1)]);                     % constant term 
    case 1
        pars = regress(xtrain(2:nobs), ...
                    [Theta, ...                                     % state term 
                    Const, ...                                      % initial condition 
                    ttrain(2:nobs)-ttrain(1), ...                   % constant term 
                    (ttrain(2:nobs).^2-ttrain(1).^2)/2]);           % linear term 
    case 2
        pars = regress(xtrain(2:nobs), ...
                    [Theta, ...                                     % state term 
                    Const, ...                                      % initial condition 
                    ttrain(2:nobs)-ttrain(1), ...                   % constant term 
                    (ttrain(2:nobs).^2-ttrain(1).^2)/2, ...         % linear term 
                    (ttrain(2:nobs).^3-ttrain(1).^3)/3 ]);          % quadratic term
    case 3
        pars = regress(xtrain(2:nobs), ...
                    [Theta, ...                                     % state term 
                    Const, ...                                      % initial condition 
                    ttrain(2:nobs)-ttrain(1), ...                   % constant term 
                    (ttrain(2:nobs).^2-ttrain(1).^2)/2, ...         % linear term 
                    (ttrain(2:nobs).^3-ttrain(1).^3)/3, ...         % quadratic term
                    (ttrain(2:nobs).^4-ttrain(1).^4)/4 ]);          % cubic term 
end

% step 3:closed-form solution to original ode
Aim = pars(1); Condim = pars(2); 

syms x(t)
switch poly_order
    case -1
        odes = diff(x) == Aim*x;                                % expression of ode 
    case 0
        Cim = pars(3);
        odes = diff(x) == Aim*x+Cim;                            % expression of ode        
    case 1
        Cim = pars(3); B1im = pars(4);
        odes = diff(x) == Aim*x+B1im*t+Cim;                     % expression of ode 
    case 2
        Cim = pars(3); B1im = pars(4); B2im = pars(5);
        odes = diff(x) == Aim*x+B1im*t+B2im*t^2+Cim;            % expression of ode 
    case 3
        Cim = pars(3); 
        B1im = pars(4); B2im = pars(5); B3im = pars(6);
        odes = diff(x) == Aim*x+B1im*t+B2im*t^2+B3im*t^3+Cim;	% expression of ode 
end

cond = x(ttrain(1)) == Condim;                                  % initial condition of ode 
xSol(t) = dsolve(odes, cond);

% step 4: fitting and forecasts of orginal series
xt_im_fits = eval(subs(xSol(t), t, [ttrain; ttest]));

% step 5: fiiting and multi-step-ahead forecasting errors 
mape_im = abs(xt_im_fits - [xtrain; xtest])./([xtrain; xtest])*100;

%% grey modelling-based ode model 
poly_order = poly_order+1;

% step 1: data preparation 
[nobs,nsta] = size(xtrain);
tgap = [1; diff(ttrain)]; 
ytrain = tgap .* ones(1,nsta) .* cumsum(xtrain);                  % cumulative sum 
Omega = (ytrain(1:nobs-1,:)+ytrain(2:nobs,:))/2;
Const = ones(nobs-1,1);

% step 2: parameter estimation
switch poly_order
    case 0      % d/dt y = a*y+c
        pars = regress(xtrain(2:nobs), ...
                    [Omega, ...                                   % state term 
                    Const]);                                      % constant term 
    case 1      % d/dt y = a*y+b1*t^1+c
        pars = regress(xtrain(2:nobs), ...
                    [Omega, ...                                   % state term 
                    Const, ...                                    % constant term 
                    (ttrain(1:nobs-1)+ttrain(2:nobs))/2 ]);       % linear term 
    case 2      % d/dt y = a*y+b2*t^2+b1*t^1+c
        pars = regress(xtrain(2:nobs), ...
                    [Omega, ...                                   % state term 
                    Const, ...                                    % constant term 
                    (ttrain(1:nobs-1)+ttrain(2:nobs))/2, ...      % linear term 
                    (ttrain(1:nobs-1).^2+ttrain(2:nobs).^2)/2 ]); % quadratic term
    case 3      % d/dt y = a*y+b3*t^3+b2*t^2+b1*t^1+c
        pars = regress(xtrain(2:nobs), ...
                    [Omega, ...                                   % state term 
                    Const, ...                                    % constant term 
                    (ttrain(1:nobs-1)+ttrain(2:nobs))/2, ...      % linear term 
                    (ttrain(1:nobs-1).^2+ttrain(2:nobs).^2)/2,... % quadratic term
                    (ttrain(1:nobs-1).^3+ttrain(2:nobs).^3)/2 ]); % cubic term 
    case 4      % d/dt y = a*y+b4*t^4+b3*t^3+b2*t^2+b1*t^1+c
        pars = regress(xtrain(2:nobs), ...
                    [Omega, ...                                   % state term 
                    Const, ...                                    % constant term 
                    (ttrain(1:nobs-1)+ttrain(2:nobs))/2, ...      % linear term 
                    (ttrain(1:nobs-1).^2+ttrain(2:nobs).^2)/2,... % quadratic term
                    (ttrain(1:nobs-1).^3+ttrain(2:nobs).^3)/2,... % cubic term 
                    (ttrain(1:nobs-1).^4+ttrain(2:nobs).^4)/2 ]); % forth term 
end

% step 3:closed-form solution to original ode
Agm = pars(1); Cgm = pars(2);

syms y(t)
switch poly_order
    case 0
        Condgm = (eye(nsta)-Agm)\Cgm; 
        odes = diff(y) == Agm*y+Cgm;                            % expression of ode        
    case 1
        B1gm = pars(3);
        Condgm = (eye(nsta)-Agm)\(Cgm+B1gm*ttrain(1));          % initial condiction
        odes = diff(y) == Agm*y+B1gm*t+Cgm;                     % expression of ode 
    case 2
        B1gm = pars(3); B2gm = pars(4);
        Condgm = (eye(nsta)-Agm)\(Cgm+B1gm*ttrain(1)+B1gm*ttrain(1)^2);
        odes = diff(y) == Agm*y+B1gm*t+B2gm*t^2+Cgm;            % expression of ode 
    case 3
        B1gm = pars(3); B2gm = pars(4); B3gm = pars(5);
        Condgm = (eye(nsta)-Agm)\(Cgm+B1gm*ttrain(1)+B1gm*ttrain(1)^2+B3gm*ttrain(1)^3);
        odes = diff(y) == Agm*y+B1gm*t+B2gm*t^2+B3gm*t^3+Cgm;	% expression of ode 
    case 4
        B1gm = pars(3); B2gm = pars(4); B3gm = pars(5); B4gm = pars(6);
        Condgm = (eye(nsta)-Agm)\(Cgm+B1gm*ttrain(1)+B1gm*ttrain(1)^2+B3gm*ttrain(1)^3+B3gm*ttrain(1)^4);
        odes = diff(y) == Agm*y+B1gm*t+B2gm*t^2+B3gm*t^3+B4gm*t^4+Cgm;	% expression of ode 
end

cond = y(ttrain(1)) == Condgm;      % initial value
ySol(t) = dsolve(odes, cond);

% step 4: fitting and forecasts of cumulative sum series
yt_im_fits = eval(subs(ySol(t), t, [ttrain; ttest]));

% step 5: fitting and forecasts of originla series
tgap = [1; diff([ttrain; ttest])];
xt_gm_fits = [yt_im_fits(1); diff(yt_im_fits)]./tgap; 

% step 6: fiiting and multi-step-ahead forecasting errors 
mape_gm = abs(xt_gm_fits - [xtrain; xtest])./([xtrain; xtest])*100;

%% Figure 7: comparison between integral matching and grey modelling 
figure('name','comparison')
plot(2004:2020,[xtrain; xtest],'+k','LineWidth',1.5)
hold on
plot(2004:2020,xt_gm_fits,'b-d','LineWidth',1.0)
plot(2004:2020,xt_im_fits,'r--o','LineWidth',1.0)
grid on
xlim([2003 2021]); xlabel('year')
ylabel('water supply (10^9 {m}^3)')
text(2012,80,' in-sample \leftarrow') 
text(2016,125,' \rightarrow out-of-sample') 
line([2015.5 2015.5], [0 150], 'linestyle', '-.', 'LineWidth',1.0);
legend({'true value','grey modelling','integral matching'},'Location','Northwest')
title('grey modeling and integral matching results')

%% 





