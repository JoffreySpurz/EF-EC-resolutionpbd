clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Resolution de l'equation de la chaleur avec 
%%% conditions aux limites de Dirichlet homogene, 
%%% par la methode des elements finis P1 et un
%%% maillage fixe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Maillage %%%%%%%%%%
% Fissure carree
pos_x = -0.5;
pos_y = 0.5;
epaisseur = 0.25;
longueur = 1;

% 1 Unit square
nodes = [-1 -1; 1 -1 ; 1 1; -1 1;pos_x pos_y ; pos_x+longueur pos_y; pos_x+longueur pos_y+epaisseur ;pos_x pos_y+epaisseur];
edges = {{1,2} {2,3} {3,4} {4,1} {5,6,'lr'} {6,7,'lr'} {7,8,'lr'} {8,5,'lr'}}; 
edges2 = {{1,2} {2,3} {3,4} {4,1}};
domain1 = Domain(nodes,edges);
%domain1 = Domain(nodes,edges2);

% mesh
dx = 0.1; % taille d'un element
mesh1 = Mesh(domain1,dx);

% Liste des sommets
P = mesh1.nodes;

% Liste des triangles
T = mesh1.triangles;

% Liste des points du bords
I = mesh1.boundary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Parametres %%%%%%%%%%
% Nombre de sommets
nP = length(P); 

% Nombre de triangles
nT = length(T);

% Condition initiale
U = sparse(nP, 1); 

% Conductivite thermique (W/m/K)
lambda_plaque = 50.2; % valeur de l'acier 
lambda_fiss = 0.024; %  valeur de l'air
% lambda_plaque =1;
% lambda_fiss =1;
% lambda_fiss = lambda_plaque;

% Masse volumique (Kg/m^3)
rho_plaque = 7850;% valeur de l'acier
rho_fiss = 1;% valeur ad hoc

% Capacite thermique massique (J/K/Kg)
cp_plaque = 444;% valeur de l'acier
cp_fiss = 1004;% valeur de l'air

% Beta = Rho*cp (K/K/m^3)
beta_plaque = rho_plaque*cp_plaque; % rho*cp dans la plaque
beta_fiss = rho_fiss*cp_fiss ; % rho*cp dans la fissure
% beta_plaque = 1;
% beta_fiss = 1;
% beta_fiss = beta_plaque;

% Plaque avec fissures
alpha = 1* mesh1.P0([ num2str(pos_x) '<=x']).*mesh1.P0(['x<=' num2str(pos_x+longueur)]).*mesh1.P0([num2str(pos_y) '<=y']).*mesh1.P0(['y<= ' num2str(pos_y+epaisseur)]); % Indicatrice : 1 dans la fissure 0 sinon
Beta = beta_plaque*(1+(beta_fiss-beta_plaque)/beta_plaque*alpha); % Vecteur Rho*cp sur chaque triangle
Lambda = lambda_plaque*(1+(lambda_fiss-lambda_plaque)/lambda_plaque*alpha); % Vecteur Conductivite thermique (lambda) sur les triangles

% Iteration
niter = 15; % Nombre d'iterations pour la resolution
dt = 0.1; % Pas en temps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Matrices d'iterations %%%%%%%%
%%% Raideur
% Paramètre : Rho*cp sur chaque triangle
LAMBDA = mesh1.P0(Lambda,2);
% Matrice de raideur
Kc = dt*mesh1.stiffness(LAMBDA);%matriceK(P, T, Lambda, dt);


%%% Masse
% Paramètre : Conductivite thermique (lambda) sur les triangles
BETA = mesh1.P0(Beta);
% Matrice de masse
M = mesh1.mass(BETA);
%mesh1.mass(Beta.*ones(nT,1));


%%% Membre de gauche
% Probleme AX(n) = MX(n-1) + F
A = M + Kc;


%%% Laser
Sv = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)
S = Sv *10^(-6); % Puissance volumique du laser utilisee
% Second membre
Fc = dt*S*mesh1.P1('x.^2+y.^2<0.2^2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Resolution %%%%%%%%%%%

%%% Avec laser
for i = 1:niter
    %%% Affichage
    mesh1.surf(U)
    caxis([0 0.04])
    pause(0.1)
    
    %%% Iteration
    U = Usolve(A, M, Fc, I, U);
end

% % On indique que le laser va etre coupe
% mesh1.surf(zeros(nP,1));
% pause(1);
% 
% %%% Sans laser
% for i = 11:niter
%     %%% Affichage
%     mesh1.surf(U)
%     caxis([0 0.04])
%     pause(0.1)
%     
%     %%% Iteration
%     U = Usolve(A, M, zeros(nP,1), I, U);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Affichage au temps final %%%%%%%%%
% mesh1.surf(U)
% caxis([0 1000])
% 
% figure;
% mesh1.plot

% figure;
% mesh1.plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%