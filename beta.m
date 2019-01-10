function [ Beta ] = beta( nT )
%%%%% Plaque %%%%%
rho_plaque = 7850;% masse volumique (Kg/m^3) : valeur de l'acier
cp_plaque = 444;% capacite thermique massique (J/K/Kg) : valeur de l'acier
beta_plaque = rho_plaque*cp_plaque; % rho*cp dans la plaque
%%%%%%%%%%%%%%%%%%


%%%%% Fissure %%%%%
rho_fiss = 100;% masse volumique (Kg/m^3) : valeur ad hoc
cp_fiss = 1004;% capacite thermique massique (J/K/Kg) : valeur de l'air
beta_fiss = rho_fiss*cp_fiss ; % rho*cp dans la fissure

% Position
ListeIndice = [ ];%35
%%%%%%%%%%%%%%%%%%


%%%%% Beta %%%%%%
% Initialisation de la plaque
Beta = beta_plaque*ones(nT, 1); 

% Prise en compte de la fissure
Beta(ListeIndice) = beta_fiss;
%%%%%%%%%%%%%%%
end

