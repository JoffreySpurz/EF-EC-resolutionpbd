function [ Flaser ] = flaser( nP, nT, T )
%%%%% Laser %%%%%
S = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)

% Position
ListeIndice = [ nT-1 ];
%%%%%%%%%%%%%%%%%%


%%%%% Flaser %%%%%%
% Initialisation
Flaser = zeros(nP,1);

%%% valeur du laser sur chaque triangle
% Initialisation de la plaque
Flaser_Vec = zeros(nT, 1); 

% Prise en compte de la fissure
Flaser_Vec(ListeIndice) = S;

%%% Iteration sur les triangles
for i=1:nT
    % On calcule Flaser(ik)
    for k=1:3
            Flaser(T(i,k)) = Flaser(T(i,k)) + Flaser_Vec(i);
    end
end
%%%%%%%%%%%%%%%
end

