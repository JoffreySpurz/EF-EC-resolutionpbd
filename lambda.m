function [ Lambda ] = lambda( nT )
%%%%% Plaque %%%%%
lambda_plaque = 50.2; % conductivite thermique (W/M/K) : valeur de l'acier 
%%%%%%%%%%%%%%%%%%


%%%%% Fissure %%%%%
lambda_fiss = 0.024; % conductivite thermique (W/M/K) : valeur de l'air

% Position
ListeIndice = [  ]; %35
%%%%%%%%%%%%%%%%%%


%%%%% Lambda %%%%%%
% Initialisation de la plaque
Lambda = lambda_plaque*ones(nT, 1); 

% Prise en compte de la fissure
Lambda(ListeIndice) = lambda_fiss;
%%%%%%%%%%%%%%%
end