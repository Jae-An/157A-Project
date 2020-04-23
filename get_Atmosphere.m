function [rho,P,T] = get_Atmosphere(z,h0)
% Atmospheric properties from NASA:
% https://www.grc.nasa.gov/www/k-12/airplane/atmos.html
% Units: (slug/ft^3), (lbf/ft^2), (F)
% h0 is launch altitude [ft]
    h = z+h0;

    if h<36152
        T = 59-(0.00356*h);
        P = 2116*(((T+459.7)/518.6)^5.256);
    end
    if h>=36152 && h<82345
        T = -70;
        P = 473.1*exp(1.73-(0.000048*h));
    end
    if h>=82345
        T = -205.05+(0.00164*h);
        P = 51.97*(((T+459.7)/389.98)^-11.388);
    end
    
    rho = P/(1718*(T+459.7));