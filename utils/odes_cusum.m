function [ Agm, Cgm, Condgm, Xt_fits ] = odes_cusum(t, Xt, tspan, solver)
% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

% Grey modelling for Cusum series
%   xt = [xt1(t_1:t_n), xt2(t_1:t_n)] 
%   tspan: time index including [ttrain; ttest]
%   model formula: d/dt yt = A*yt + c

[nobs, nsta] = size(Xt);                    % size of Xt
%% cumulative sum
h = [1; diff(t)]; 
Yt =  cumsum([h h].*Xt);

%% parameter estimation: two column with two regression 
Const = ones(nobs-1,1);
Theta = (Yt(1:nobs-1,:)+Yt(2:nobs,:))/2;

Xi = zeros(nsta, nsta+1);
for ind_sta = 1:nsta
    pars = regress(Xt(2:nobs,ind_sta), [Theta, Const]);
    Xi(ind_sta,:) = pars';
end

Agm = Xi(:,1:nsta); Cgm = Xi(:,1+nsta);
Condgm = (eye(nsta)-Agm)\Cgm;

if solver == 1
    %% method II: the expression of the close-form solution (FAST)
    % obtaiin the closed-form solutions 
    % syms y1(t) y2(t)
    % syms a11 a12 a21 a22 c1 c2 cond1 cond2
    % y = [y1; y2];
    % Agm = [a11 a12; a21 a22]; Cgm = [c1; c2]; Condgm = [cond1; cond2];
    % odes = diff(y) == Agm*y+Cgm; 
    % cond = y(0) == Condgm;
    % [y1Sol(t), y2Sol(t)] = dsolve(odes, cond);
    % y1Sol(t)
    % y2Sol(t)

    % % step 2: fitting and forecasting values of Cusum and Orgin series
    a11 = Agm(1,1); a12 = Agm(1,2); a21 = Agm(2,1); a22 = Agm(2,2); 
    c1 = Cgm(1); c2 = Cgm(2); cond1 = Condgm(1); cond2 = Condgm(2);

    nobs = length(tspan);
    Yt_fits = zeros(nobs, nsta);
    for iter_t = 1:nobs
        t = tspan(iter_t);
        Yt_fits(iter_t,:) = ...
        [exp((t*(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*(a22/a21 - (a11/2 + a22/2 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)/2)/a21)*((4*a21*c1 - 2*a11*c2 + 2*a22*c2 - 2*c2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) - a11^2*cond2 + a22^2*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21) + 2*a11*a21*cond1 + 2*a21*a22*cond1 - 2*a21*cond1*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) - 2*a22*cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/(2*(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)) - (exp((t*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/2 - (a22*t)/2 - (a11*t)/2)*(2*c2 - (2*a11*c2 - 4*a21*c1 - 2*a22*c2 + 2*c2*(a11 + a22))/(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))) - exp((t*(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*((4*a21*c1 - 2*a11*c2 + 2*a22*c2 + 2*c2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) - a11^2*cond2 + a22^2*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21) + 2*a11*a21*cond1 + 2*a21*a22*cond1 + 2*a21*cond1*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) + 2*a22*cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/(2*(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)) - (exp(- (a11*t)/2 - (a22*t)/2 - (t*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/2)*(2*c2 - (2*a11*c2 - 4*a21*c1 - 2*a22*c2 + 2*c2*(a11 + a22))/(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))*(a22/a21 - (a11/2 + a22/2 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)/2)/a21), ...
         exp((t*(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*((4*a21*c1 - 2*a11*c2 + 2*a22*c2 + 2*c2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) - a11^2*cond2 + a22^2*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21) + 2*a11*a21*cond1 + 2*a21*a22*cond1 + 2*a21*cond1*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) + 2*a22*cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/(2*(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)) - (exp(- (a11*t)/2 - (a22*t)/2 - (t*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/2)*(2*c2 - (2*a11*c2 - 4*a21*c1 - 2*a22*c2 + 2*c2*(a11 + a22))/(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))) - exp((t*(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*((4*a21*c1 - 2*a11*c2 + 2*a22*c2 - 2*c2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) - a11^2*cond2 + a22^2*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21) + 2*a11*a21*cond1 + 2*a21*a22*cond1 - 2*a21*cond1*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2) - 2*a22*cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/(2*(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)) - (exp((t*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))/2 - (a22*t)/2 - (a11*t)/2)*(2*c2 - (2*a11*c2 - 4*a21*c1 - 2*a22*c2 + 2*c2*(a11 + a22))/(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2))) ...
        ];
    end
else
    %% method I: close-form soluion to ODEs: d/dt y = A*y + c (SLOW)
    % step 1: dsolve to obtain the closed-from solution 
    syms y1(t) y2(t)
    y = [y1; y2];
    odes = diff(y) == Agm*y+Cgm; 
    cond = y(0) == Condgm; 
    [y1Sol(t), y2Sol(t)] = dsolve(odes, cond);

    % step 2: fitting and forecasting values of Cusum and Orgin series
    Yt_fits = [eval(subs(y1Sol(t), t, tspan)), ...
                eval(subs(y2Sol(t), t, tspan))];
end
       
%% inverse cumulative sum 
hh = [1; diff(tspan)]; 
[nobs, nsta] = size(Yt_fits);
Xt_fits = zeros(nobs, nsta);
for ind_col=1:nsta
    Xt_fits(:,ind_col) = [Condgm(ind_col); ...
        Yt_fits(2:nobs,ind_col)-Yt_fits(1:nobs-1,ind_col)]./hh;        
end

end

