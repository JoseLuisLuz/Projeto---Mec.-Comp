
%------------------------------------------------------------------------------------------------------------------
% - Projeto Mec�nica Computacional - Licenciatura Engenharia Mec�nica - 3� ano
% - Lemec 21/22
% - An�lise Linear de Tens�o plana numa placa com entalhe;
% - Trabalho realizado por:
%                           - Jos� Luz - 103489;
%                           - Miguel Cola�o - 103370;
%                           - Miguel Vieira - 103359;
% - Orientador: Professor Pedro Areias
%
% - Fun��o das regras de integra��p num�rica para triangulos. Este
% c�digo foi retirado dos documentos proporcionados pelo professor Leonel
% Fernandes
% 
%------------------------------------------------------------------------------------------------------------------


function [xp wp]=GenipT (nip)
% pesquisar: S. Deng quadrature formulas in two dimensions matlab
if (nip == 3)   % regra de 3 pts e grau 2
xp=[0.5 0; 0.5 0.5;0 0.5];   % (3x2)
wp=[1 ; 1; 1]/6;
end
if (nip == 4) % regra de 4 pts e grau 3
xp=[1/3 1/3; 0.2 0.2;0.6 0.2;0.2 0.6]; % (4x2)
wp=[-27 ; 25; 25; 25]/96;
end
if (nip == 6) % regra de 6 pts e grau 4
xp=[0.44594849091597 0.44594849091597 ; % (6x2)
0.44594849091597 0.10810301816807 ;
0.10810301816807 0.44594849091597 ;
0.09157621350977 0.09157621350977 ;
0.09157621350977 0.81684757298046 ;
0.81684757298046 0.09157621350977 ];

wp=[ 0.22338158967801; 0.22338158967801; 0.22338158967801;
 0.10995174365532 ; 0.10995174365532; 0.10995174365532]/2;
end
if (nip == 7) % regra de 7 pts e grau 5
xp=[0.33333333333333 0.33333333333333 ;  % (7x2)
    0.47014206410511 0.47014206410511 ;
    0.47014206410511 0.05971587178977 ;
    0.05971587178977 0.47014206410511 ;
    0.10128650732346 0.10128650732346 ;
    0.10128650732346 0.79742698535309 ; 
    0.79742698535309 0.10128650732346];
wp=[ 0.22500000000000; 0.13239415278851;
 0.13239415278851; 0.13239415278851;
 0.12593918054483; 0.12593918054483;
 0.12593918054483]/2;
end
end