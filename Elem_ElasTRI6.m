
%------------------------------------------------------------------------------------------------------------------
% - Projeto Mecânica Computacional - Licenciatura Engenharia Mecânica - 3º ano
% - Lemec 21/22
% - Análise Linear de Tensão plana numa placa com entalhe;
% - Trabalho realizado por:
%                           - José Luz - 103489;
%                           - Miguel Colaço - 103370;
%                           - Miguel Vieira - 103359;
% - Orientador: Professor Pedro Areias
%
% - Função do cálculo da matriz rigidez e vetor de forças para um triangulo T6. Este
% código foi retirado dos documentos proporcionados pelo professor Leonel
% Fernandes
% 
%------------------------------------------------------------------------------------------------------------------


function [Ke fe]=Elem_ElasTRI6 (XN,C,he,fL)
%   Matriz XN(6,2) contem as coordenadas locais deste TRI de 6 nos
%   fL - intensidade da força uniforme  
%   inicializar Ke
 Ke = zeros(12,12);
 A = zeros(6,2);	% matriz auxiliar A(6x2) das derivadas parciais das
%	 funcoes de forma em (x,y)
 F = zeros(6,1);	% vector auxiliar F(6,1) das funcoes de forma
%--------------------------------------------------------------------------
%   gerar os pontos de integracao
nip = 7;
[xp wp]=GenipT (nip);

%   percorrer os pontos de integracao
for ip=1:nip;	%	ciclo para os pontos de integracao

csi = xp(ip,1);
eta = xp(ip,2);
%
%   para cada ponto dado, calcular:   
%   a matriz A (6x2) das derivadas parciais das funcoes de forma em (x,y);  
%   funcoes de forma do TRI-6, no vector F(6x1);
%   jacobiano da transformacao, Detj.
%-------------------------------------------------------
[A F Detj]=Shape_N_Der6 (XN,csi,eta);
%--------------------------------------------------------------------------
%      forma a matriz B (3x12) para elasticidade plana
B=[A(1,1) 0 A(2,1) 0 A(3,1) 0 A(4,1) 0 A(5,1) 0 A(6,1) 0;0 A(1,2) 0 A(2,2)...
   0 A(3,2) 0 A(4,2) 0 A(5,2) 0 A(6,2); A(1,2) A(1,1) A(2,2) A(2,1) A(3,2)...
   A(3,1) A(4,2) A(4,1) A(5,2) A(5,1) A(6,2) A(6,1)];
%--------------------------------------------------------------------------
%      forma a matriz Psi (2x12) para elasticidade plana
Psi=[F(1) 0 F(2) 0 F(3) 0 F(4) 0 F(5) 0 F(6) 0;0 F(1) 0 F(2) 0 F(3) 0 ...
     F(4) 0 F(5) 0 F(6)];
%--------------------------------------------------------------------------
%   7) peso transformado
wip = he*wp(ip)*Detj;
%   8) calcular produto B'*C*B, pesar e somar a Ke
Ke = Ke + wip*B'*C*B ;
%   9) %   vector de cargas: somar fe aqui quando tiver de ser ;)
%wipf = fL*wip ;
%fe = fe + wipf*Psi ;
fe = zeros(12,1);    %   nao ha body forces %%???? Juro que esta parte me deixa perplexo
%	proximo ponto
end	%   fim de ciclo de integracao
end


function [xp wp]=GenipT (nip)
% pesquisar: S. Deng quadrature formulas in two dimensions matlab
if (nip == 3)   % regra de 3 pts e grau 2
xp=[0.5 0; 0.5 0.5;0 0.5];   % (3x2)
wp=[1 ; 1; 1]/6;
elseif (nip == 4) % regra de 4 pts e grau 3
xp=[1/3 1/3; 0.2 0.2;0.6 0.2;0.2 0.6]; % (4x2)
wp=[-27 ; 25; 25; 25]/96;
elseif (nip == 6) % regra de 6 pts e grau 4
xp=[0.44594849091597 0.44594849091597 ; % (6x2)
0.44594849091597 0.10810301816807 ;
0.10810301816807 0.44594849091597 ;
0.09157621350977 0.09157621350977 ;
0.09157621350977 0.81684757298046 ;
0.81684757298046 0.09157621350977 ];

wp=[ 0.22338158967801; 0.22338158967801; 0.22338158967801;
 0.10995174365532 ; 0.10995174365532; 0.10995174365532]/2;
elseif (nip == 7) % regra de 7 pts e grau 5
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

function [B psi Detj]=Shape_N_Der6 (XN,csi,eta)
%----------------------------------------------------------------
%   Matriz XN(6,2) contem as coordenadas locais deste triangulo de 6 nos
%----------------------------------------------------------------
psi=zeros(6,1) ;
%   para cada ponto dado, calcular
%   1) funcoes de forma do tri-6, vector psi (6x1)
psi(1) = (1-csi-eta)*(1-2*csi-2*eta);
psi(2) = csi*(2*csi-1);
psi(3) = eta*(2*eta-1);
psi(4) = 4*(1-csi-eta)*csi;
psi(5) = 4*csi*eta;
psi(6) = 4*(1-csi-eta)*eta;

%   2) derivadas parciais em (csi,eta), Matriz Dpsi(6x2)
Dpsi(1,1) = 4*csi+4*eta-3;
Dpsi(2,1) = 4*csi-1;
Dpsi(3,1) = 0;
Dpsi(4,1) = 4 -8*csi -4*eta;
Dpsi(5,1) = 4*eta;
Dpsi(6,1) = -4*eta;
%
Dpsi(1,2) = 4*csi+4*eta-3;
Dpsi(2,2) = 0;
Dpsi(3,2) = 4*eta-1;
Dpsi(4,2) = -4*csi;
Dpsi(5,2) = 4*csi;
Dpsi(6,2) = 4 -4*csi -8*eta;
%   3) derivadas parciais da matriz jacobiana (2x2) de x e y
jaco = XN'*Dpsi ;
%   4) jacobiano da transformacao
Detj = det(jaco) ;
%   5) derivadas parciais da transformacao inversa (csi,eta) em funcao de
%     (x,y), matriz (2x2)
Invj = inv(jaco) ;
%   6) finalmente a matriz B (6x2) das derivadas parciais das funcoes de 
%      forma em (x,y)
B = Dpsi*Invj ;
%----------------------------------------------------------------
end