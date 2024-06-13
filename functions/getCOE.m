function [a, e, i, Omega, w, theta] = getCOE(r_vect, v_vect, mu)
tol = 1e-10;
% 1) Calculate distance of s/c from foci
r = norm(r_vect);
% 2) Calculate velocity magnitude
v = norm(v_vect);
% 3) Calculate specific angular momentum
h_vect = cross(r_vect, v_vect);
h = norm(h_vect);
% 4) Calculate Eccentricty vector and eccentricty
e_vect = 1/mu*(cross(v_vect, h_vect) - mu/r*r_vect);
e = norm(e_vect);
% 5) Calculate semi major axis
a = mu/(2*mu/r - v^2);
% 6) Calculate Inclination
i = acos(h_vect(3)/h);
% 7) Calculate Nodal vector
N_vect = cross([0;0;1], h_vect);
% 8) Calculate magnitude of nodal vector
N = norm(N_vect);
% 9) Calculate Right ascension of ascending node
if N~= 0
    Omega = acos(N_vect(1)/N);
    if N_vect(2) < 0
        Omega = 2*pi - Omega;
    end
else
    Omega = 0;
end
% 9) Calculate argument of perigee
if N~= 0
    if e> tol
        w = acos(dot(N_vect,e_vect)/(N*e));
        if e_vect(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end
% 10) Calculate true anomaly
if e > tol
    theta = acos(dot(e_vect, r_vect)/(e*r));
    if dot(r_vect, v_vect) < 0
        theta = 2*pi - theta;
    end
else
    cp = cross(N_vect, r_vect);
    theta = acos(dot(N_vect, r_vect)/(r*N));
    if cp(3) < 0
        theta = 2*pi - theta;
    end
end
end