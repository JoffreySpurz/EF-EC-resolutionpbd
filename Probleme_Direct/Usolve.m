%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B : Matrice (nP x nP) de la forme A\Mc
% Fc : Vecteur (nP x 1) du laser
% U : matrice de la temperature dans la plaque
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ U ] = Usolve(B, Fc, nP, niter, niter_laser)
%%% Condition initiale %%%
U = zeros(nP, niter); 
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Si niter_laser est defini
if nargin == 5
    
%%% Resolution %%%
    % Avec laser
    for i=2:niter_laser
        U(:,i) = B*(U(:,i-1)+Fc);
    end
    
    % Sans laser
    for i=niter_laser+1:niter
        U(:,i) = B*U(:,i-1);
    end
%%%%%%%%%%%%%%%%%%%%%


% Sinon
else 
    
%%% Resolution %%%
    % Avec laser
    for i=2:niter
        U(:,i) = B*(U(:,i-1)+Fc);
    end
%%%%%%%%%%%%%%%%%%%%%
end

end

