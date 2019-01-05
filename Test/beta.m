function [ Beta ] = beta( nP, nT, T )
%%%%% Plaque %%%%%
rho_plaque = 7850;% masse volumique (Kg/m^3) : valeur de l'acier
cp_plaque = 444;% capacite thermique massique (J/K/Kg) : valeur de l'acier
beta_plaque = rho_plaque*cp_plaque; % rho*cp dans la plaque
%%%%%%%%%%%%%%%%%%


%%%%% Fissure %%%%%
rho_fiss = 1000;% masse volumique (Kg/m^3) : valeur ad hoc
cp_fiss = 1004;% capacite thermique massique (J/K/Kg) : valeur de l'air
beta_fiss = rho_fiss*cp_fiss ; % rho*cp dans la fissure

% Position
ListeIndice = [ 35 ];
%%%%%%%%%%%%%%%%%%


%%%%% Beta %%%%%%
% Initialisation de la matrice
Beta = zeros(nP,nP);

%%% valeur de beta sur chaque triangle
% Initialisation de la plaque
Beta_Vec = beta_plaque*ones(nT, 1); 

% Prise en compte de la fissure
Beta_Vec(ListeIndice) = beta_fiss;

%%% Iteration sur les triangles
for i=1:nT
    % On calcule Beta(ik,il)
    for k=1:3
        for l=1:3 
            Beta(T(i,k), T(i,l)) = Beta(T(i,k), T(i,l)) + Beta_Vec(i);
        end
    end
end
%%%%%%%%%%%%%%%
end

