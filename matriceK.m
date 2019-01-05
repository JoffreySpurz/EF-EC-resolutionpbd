%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P : liste des sommets du maillage, pour la ieme ligne : Pi = [Pi1 Pi2]
% T : liste des triangles du maillage, pour la ieme ligne : Ti = [nPi1 nPi2
% nPi3]
% Lambda : vecteur des lambda sur chaque triangle
% dt : pas en espace
% K : la matrice de rigidite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ K ] = matriceK( P, T, Lambda, dt )
% Nombre de triangles
NT = length(T);

% Nombre de sommets
NP = length(P);

% Initialisation creuse de K
K = sparse(NP, NP);

% Iteration sur les triangles
for i=1:NT
    % On recupere Pi
    P1 = P(T(i,1),:);
    P2 = P(T(i,2),:);
    P3 = P(T(i,3),:);
    
    % On definit dn
    [S,gn] = geo(P1', P2', P3');
        
    % On recupere a, b
    [a, b] = base(P1', P2', P3');
    
    % On calcule Aik,il
    for k=1:3
        for l=1:3
            K(T(i,k), T(i,l)) = K(T(i,k), T(i,l)) + S*(a(1:2, k)' * a(1:2, l)) * Lambda(i) * dt;
        end
    end
end

end

