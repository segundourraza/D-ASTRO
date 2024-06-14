function [r_vect, v_vect] = getECIframe(R, lon, lat, V, gamma, A)

% Orbital position
r_x = cos(lat).*cos(lon).*R;
r_y = cos(lat).*sin(lon).*R;
r_z = sin(lat).*R;
r_vect = [r_x, r_y, r_z];

% Velocity vector
v_x = V.*sin(gamma);
v_y = V.*cos(gamma).*sin(A);
v_z = V.*cos(gamma).*cos(A);
v_vect = [v_x, v_y, v_z];
