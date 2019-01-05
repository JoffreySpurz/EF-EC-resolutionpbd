function [ Lambda ] = lambda( nP, nT, T )
%%%%% Plaque %%%%%
lambda_plaque = 50.2; % conductivite thermique (W/M/K) : valeur de l'acier 
%%%%%%%%%%%%%%%%%%


%%%%% Fissure %%%%%
lambda_fiss = 0.024; % conductivite thermique (W/M/K) : valeur de l'air

% Position
ListeIndice = [ 35 ];
%%%%%%%%%%%%%%%%%%


%%%%% Lambda %%%%%%
% Initialisation de la matrice
Lambda = zeros(nP,nP);

%%% valeur de lambda sur chaque triangle
% Initialisation de la plaque
Lambda_Vec = lambda_plaque*ones(nT, 1); 

% Prise en compte de la fissure
Lambda_Vec(ListeIndice) = lambda_fiss;

%%% Iteration sur les triangles
for i=1:nT
    % On calcule Lambda(ik,il)
    for k=1:3
        for l=1:3 
            Lambda(T(i,k), T(i,l)) = Lambda(T(i,k), T(i,l)) + Lambda_Vec(i);
        end
    end
end
%%%%%%%%%%%%%%%
end