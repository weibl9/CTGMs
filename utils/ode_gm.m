function dx = ode_gm(t, x, beta)
% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

    dx = beta(1)*x + beta(2)*t^2 + beta(3)*t + beta(4);
    
end