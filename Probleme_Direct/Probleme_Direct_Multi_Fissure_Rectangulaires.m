clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Resolution de l'equation de la chaleur avec 
%%% conditions aux limites de Dirichlet homogene, 
%%% par la methode des elements finis P1 et un
%%% maillage fixe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Maillage %%%%%%%%%%
% Fissure carree
nb_fiss = 2;% Nombre de fissures dans la plaque
Pos_x = [-0.5 0.35];% Liste des abscisses
Pos_y = [0.25 -0.75];% Liste des ordonnees
Longueur_y = [0.1 0.75];% Liste des longueurs selon y
Longueur_x = [1 0.1];% Liste des longueurs selon x

% 1 Unit square
nodes = [-1 -1; 1 -1 ; 1 1; -1 1];
edges = {{1,2} {2,3} {3,4} {4,1}};

for i=1:nb_fiss
    % Parametres de la fissure i
    pos_x = Pos_x(i);
    pos_y = Pos_y(i);
    longueur_x = Longueur_x(i);
    longueur_y = Longueur_y(i);
    
    % Mise a jour du maillage
    nodes = [nodes ; pos_x pos_y ; pos_x+longueur_x pos_y; pos_x+longueur_x pos_y+longueur_y ;pos_x pos_y+longueur_y];
    index = length(edges); % On recupere la longueur de la cellule
    edges{1,index+1} = {index+1,index+2,'lr'};  
    edges{1,index+2} = {index+2,index+3,'lr'};
    edges{1,index+3} = {index+3,index+4,'lr'};
    edges{1,index+4} = {index+4,index+1,'lr'};
end

domain1 = Domain(nodes,edges);


% mesh
dx = 0.05; % taille d'un element
mesh1 = Mesh(domain1,dx,'save');

% Liste des sommets
P = mesh1.nodes;

% Liste des triangles
T = mesh1.triangles;

% Liste des points du bords
I = mesh1.boundary([1 2 3 4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Parametres %%%%%%%%%%
% Nombre de sommets
nP = length(P); 

% Nombre de triangles
nT = length(T);

% Conductivite thermique (W/m/K)
lambda_plaque = 50.2; % valeur de l'acier 
lambda_fiss = 0.024; %  valeur de l'air

% Masse volumique (Kg/m^3)
rho_plaque = 7850;% valeur de l'acier
rho_fiss = 1;% valeur ad hoc

% Capacite thermique massique (J/K/Kg)
cp_plaque = 444;% valeur de l'acier
cp_fiss = 1004;% valeur de l'air

% Beta = Rho*cp (K/K/m^3)
beta_plaque = rho_plaque*cp_plaque; % rho*cp dans la plaque
beta_fiss = rho_fiss*cp_fiss ; % rho*cp dans la fissure

% Plaque avec fissures
alpha = mesh1.P0(0);

for i=1:nb_fiss
    % Parametres de la fissure i
    pos_x = Pos_x(i);
    pos_y = Pos_y(i);
    longueur_x = Longueur_x(i);
    longueur_y = Longueur_y(i);
    
    % Indicatrice fissure (1 dans la fissure 0 sinon)
    alpha = alpha + mesh1.P0([ num2str(pos_x) '<=x']).*mesh1.P0(['x<=' num2str(pos_x+longueur_x)]).*mesh1.P0([num2str(pos_y) '<=y']).*mesh1.P0(['y<= ' num2str(pos_y+longueur_y)]);
end
alpha = 1*alpha;

Beta = (beta_plaque/lambda_plaque)*(1+(beta_fiss-beta_plaque)/beta_plaque*alpha); % Vecteur Rho*cp sur chaque triangle
Lambda = 1+(lambda_fiss/lambda_plaque-1)*alpha; % Vecteur Conductivite thermique (lambda) sur les triangles

% Iteration
niter = 1000; % Nombre d'iterations pour la resolution
dt =10^(-1); % Pas en temps
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
% Probleme AX(n) = M( X(n-1) + F )
A = M + Kc;

%%% Prise en compte des CL %%%
% Nombre de points du bord
nI = length(I);

% CL de Dirichlet
A(I,:) = 0;
A(I,I) = speye(nI);
M(I,:) = 0;
M(I,I) = speye(nI);

% Probleme X(n) = A\M ( X(n-1) + F ) = B ( X(n-1) + F )
B = A\M;


%%% Laser
Sv = 8*10^9 ; % Puissance volumique du laser applique a la plaque (W/m^3)
S = Sv *10^(-9)/lambda_plaque; % Puissance volumique du laser utilisee
% Second membre
Fc = dt*S*mesh1.P1('exp(-(x.^2+y.^2)/0.05)');% Gaussien
%Fc=dt*S*mesh1.P1('x.^2+y.^2<0.2^2');% Cercle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Resolution %%%%%%%%%%%
%%% Avec laser
U = Usolve(B, Fc, nP, niter);

%%% Save
save('Probleme_direct_multi_fissure_rectangulaires','U', 'mesh1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%