%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pi est un sommet du triangle
% S surface du triangle
% G barycentre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ S, G ] = geo( P1, P2, P3 )
% On définit dn 
dn = abs(det( [P2-P1 P3-P1] ));

% Surface d'un triangle
S = dn/2;

% Barycentre d'un triangle
G = (1/3)*(P3+ P2+ P1);
end

