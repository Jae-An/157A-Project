function [Cd_f] = get_Cdf(rocket, v, M, nu, L, K)
    
    %% Geo
    % Body
    d = rocket.geo.body.D; % max body diameter
    L_nose = rocket.geo.nose.L;
    S_b = pi*d*(0.5*rssq([d/2,L_nose]) + L);    % Body

    % Fins
    t = rocket.geo.fins.t;
    C_r = rocket.geo.fins.RC;
    h_t = 0.5; % conservative over square fin???
    C_t = rocket.geo.fins.TC;
    N = rocket.geo.fins.N;
    S_f = rocket.geo.fins.SS * (rocket.geo.fins.RC + rocket.geo.fins.TC);

    %% Calcs
    Cf_final = zeros(2,1); % only (nose+)body, fins
    L_char = [L, C_r]; % characteristic lengths
    for i = 1:2
        Rn_str = (v*L_char(i)/(12*nu))...
               * (1 + 0.0283*M - 0.043*M^2 + 0.2107*M^3 - 0.03829*M^4 + 0.002709*M^5);
        Cf_str = 0.037036 * Rn_str^-0.155079;
        Cf = Cf_str * (1 + 0.00798*M - 0.1813*M^2 + 0.0632*M^3 - 0.00933*M^4 + 0.000549*M^5);
        Cf_str_term = (1.89 + 1.62*log10(L/K))^-2.5;
        Cf_term = Cf_str_term / (1 + 0.2044*M^2);
        Cf_final(i) = max([Cf, Cf_term]);
    end
    
    % Body
    Cd_f_body = Cf_final(1)*(1 + 60/(L/d)^3 + 0.0025*(L/d))*(4*S_b/(pi*d^2));
    
    % Fins
    Rn = v*C_r / (12*nu);
    lambda = C_t / C_r;
    Cf_lambda = Cf_final(2) * (log10(Rn)^2.6 / (lambda^2 - 1))...
              * ((lambda^2 / log10(Rn*lambda)^2.6) - (1 / log10(Rn)^2.6)...
                 + 0.5646*((lambda^2 / log10(Rn*lambda)^3.6) - (1 / log10(Rn)^3.6)));
    Cd_f_fins = Cf_lambda * (1 + 60*(t/C_r)^4 + 0.8*(1+5*h_t^2)*(t/C_r)) * 4*N*S_f/(pi*d^2);         

        
    % Total (need to add protuberance, excrescency drag)
    Cd_f = Cd_f_body + 1.04*Cd_f_fins;

end