
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
% - Função do cálculo da matriz rigidez e vetor de forças para um triangulo T3. Este
% código foi retirado dos documentos proporcionados pelo professor Leonel
% Fernandes
% 
%------------------------------------------------------------------------------------------------------------------


function [Ke fe]=Elem_CST (x1,y1,x2,y2,x3,y3,C,he,el)
%--------------------------------------------------------------------------
%   C  - Matriz Elastica para Elasticidade Plana
%	el - element load - carga uniforme - variavel de entrada
%--------------------------------------------------------------------------
%       calcula o dobro da area sinalizada do triangulo
Ae2 = (x2 -x1)*(y3 -y1) -(y2 -y1)*(x3 -x1);

if Ae2 < 0
    %disp('els com area negativa a ser corrigida')
    %Ae2 = -1*Ae2;
    Ae2 = (x2 -x3)*(y1 -y3) -(y2 -y3)*(x1 -x3);
end

%       calculo da area do triangulo
Ae= Ae2/2;	%	area sinalizada
%       derivadas parciais das funcoes de forma
d1dx = (y2-y3)/Ae2;
d1dy = (x3-x2)/Ae2;
d2dx = (y3-y1)/Ae2;
d2dy = (x1-x3)/Ae2;
d3dx = (y1-y2)/Ae2;
d3dy = (x2-x1)/Ae2;
%--------------------------------------------------------------------------
%       forma a matriz B (3x6) para elasticidade plana
B=zeros(3,6);
%       matriz B (3x6) para elasticidade plana
B=[d1dx 0 d2dx 0 d3dx 0;0 d1dy 0 d2dy 0 d3dy;d1dy d1dx d2dy d2dx d3dy d3dx];
Ke = he*Ae*B'*C*B;
%   vector de cargas 
fe(1:6,1)= 0;    %   nao ha body forces
end