function [ Flaser ] = flaser( nT )
%%%%% Laser %%%%%
S = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)

% Position
ListeIndice = [ nT-1 ];
%%%%%%%%%%%%%%%%%%


%%%%% Flaser %%%%%%
% Initialisation
Flaser = zeros(nT,1);

% Prise en compte de la fissure
Flaser(ListeIndice) = S;
%%%%%%%%%%%%%%%
end

