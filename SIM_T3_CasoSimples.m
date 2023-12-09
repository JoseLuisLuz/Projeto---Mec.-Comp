
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
% - Simulação do caso simples (quadrado sem entalhe 10mmx10mm) através de MEF por Triângulo de 3 Nós
% - Partes do código foram retiradas do código proporconado nas aulas do
% professor Leonel Fernandes
%------------------------------------------------------------------------------------------------------------------


%------------------------------------------------------------------------------------------------------------------
                %% Comandos de Limpeza
%------------------------------------------------------------------------------------------------------------------

clear all
close all

%------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------------
                %% Especificações do Elemento
%------------------------------------------------------------------------------------------------------------------


%  Dimensões
h = 10;
esp = 1;

% Propriedades do material
E = 70*10^3;
niu = 0.3;

% Relação constitutiva

cons=E/(1 - niu*niu);
C=cons*[1 niu 0;niu 1 0;0 0 (1-niu)/2];


% Coordenadas e nós dos elementos no caso simples

x = [-5; 5; 5; -5];
y = [-5; -5; 5; 5];

conects= [1 1 2 3; 2 3 4 1];
coords= [ [1:4]', x, y];

x=coords(:,2);
y= coords(:,3);

%------------------------------------------------------------------------------------------------------------------
                %% Representação gráfica do elemento 
%------------------------------------------------------------------------------------------------------------------


FIG1= figure (1);
Nelt=size(conects,1);

%------------------------------------------------------------------------------------------------------------------

for i=1: Nelt

    no1 = conects(i,2);
    no2 = conects(i,3);
    no3 = conects(i,4); 


    edofs = [no1 no2 no3 no1];
    fill (x(edofs),y(edofs),'green'); hold on
    plot(x(edofs),y(edofs),'black'); hold on

    % Possibilidade de representar nos centros o número com o elemento

     % cenx= (x(no1) + x(no2) + x(no3))/3;
     % ceny= (y(no1) + y(no2) + y(no3))/3;
     % text(cenx, ceny, num2str(i)); hold on

end

%------------------------------------------------------------------------------------------------------------------

% plot da geometria

plot(x,y,'.'); hold on
axis([-10 10 -10 10]);
title('Representação Inicial da Malha')

% Possibilidade de representar nos nós o seu número

%text(x, y, num2str(coords(:,1))); hold off

hold off

%------------------------------------------------------------------------------------------------------------------
                %% Assemblagem dos elementos na matriz Global
%------------------------------------------------------------------------------------------------------------------

Nnds = size(coords,1);        % número de nós
Neqs = 2*Nnds;            % número total de equacões ou graus de liberdade

% Definição inicial da matriz global

Kg = zeros(Neqs,Neqs);     % Matriz Kg
fg = zeros(Neqs,1);	% Vetor de forças volúmicas

%------------------------------------------------------------------------------------------------------------------

% Ciclo para os definição da matriz local de cada elemento através da
% função Elem_CST e assemblagem da matriz global

for i=1:Nelt

    no1 = conects(i,2);
    no2 = conects(i,3);
    no3 = conects(i,4);

    % Endereços na matriz global dos graus de liberdade dos 3 nós locais de
    % cada triângulo 

    idof1 = (no1-1)*2+1;
    idof2 = (no2-1)*2+1;
    idof3 = (no3-1)*2+1;

    % Vetor com os graus de liberdade do nós do elemento

    edofs = [idof1 idof1+1 idof2 idof2+1 idof3 idof3+1];

    % Matriz local do elemento e vetor força (A função foi retirada do
    % código Matlab do Professor Leonel Fernandes)

    [Ke fe] = Elem_CST (coords(no1,2), coords(no1,3),coords(no2,2),coords(no2,3),coords(no3,2),coords(no3,3),C,esp,0);   % <- carregamento nulo aqui

    % Assemblagem do elemento associando a matriz local na matriz global e
    % o vetor forças local no global

    Kg(edofs,edofs) = Kg(edofs,edofs) + Ke;
    fg(edofs,1) = fg(edofs,1) + fe;

end

%------------------------------------------------------------------------------

% Arquivo da matriz   

ka = Kg;

%------------------------------------------------------------------------------------------------------------------
                %% Definição das condições de fronteira Naturais do elemento
%------------------------------------------------------------------------------------------------------------------

% Tracões aplicadas nos nós

% Exemplo de trações por nó, vertical e horizontal respetivamente

% no=1
%     idof = (no-1)*2+2	        % Endereço do segundo grau de liberdade deste nó
%     fg(idof,1)= fg(idof,1) + 10;      % Força Vetical = 10N 

% no=2
%     idof = (no-1)*2+1	        % Endereço do primeiro grau de liberdade deste nó
%     fg(idof,1)= fg(idof,1) + 10;      % Força Horizontal = 10N 

%------------------------------------------------------------------------------------------------------------------

% Implementação para uma carga distribuida numa face da geometria

% Nós com condições naturais aplicadas

noscondat = [];

for i = 1:Nnds

    if y(i) == 5       % Aplicação de carga nos nós em y = 5   

        noscondat(end+1,:) = coords(i,1);

    end
end

Nnoscondat = size(noscondat, 1);

% Distribuição das forças pelos nós

qx = 0;
qy = 2;
fpnx = qx*h/(Nnoscondat);
fpny = qy*h/(Nnoscondat);

for i= 1: size(noscondat, 1)

        no = noscondat(i);

        iydof = (no-1)*2+2;	 %   endereco do primeiro grau de liberdade deste no
        fg(iydof,1)= fg(iydof,1) + fpny;  %Força vetical = 2N/mm (10N em cada nó)
end

for i= 1: size(noscondat, 1)

        no = noscondat(i);

        ixdof = (no-1)*2+1;	 %   endereco do primeiro grau de liberdade deste no
        fg(ixdof,1)= fg(ixdof,1) + fpnx;  %Força vetical = 2N/mm (10N em cada nó)
end


%------------------------------------------------------------------------------------------------------------------

% Arquivo da matriz global com condições naturais aplicadas para cálculo do
% equilibrio estático e reações nos apoios

Kr = Kg;
fr = fg;

ka = Kr;

%------------------------------------------------------------------------------------------------------------------
                %% Definição das condições de fronteira essenciais do elemento
%------------------------------------------------------------------------------------------------------------------

boom = 1.0e+14;

% Exemplo de condições por nó, vertical e horizontal respetivamente

% Restringir graus de liberdade em x

% no=1
%     idof = (no-1)*2+1;        % Endereço do primeiro grau de liberdade deste nó
%     Kr(idof,idof) = boom;
%     fr(idof,1)= boom*0;

% Restringir graus de liberdade y (verticais)

% no=1
%     idof = (no-1)*2+2;        % Endereço do segundo grau de liberdade deste nó
%     Kr(idof,idof) = boom;
%     fr(idof,1)= boom*0;

%------------------------------------------------------------------------------------------------------------------

% Implementação da restrição dos deslocamentos distribuida num face da geometria

% Nós com condições essenciais aplicadas

noscondess = [];

for i = 1:Nnds

    if y(i) == -5       % Aplicação da restrição dos deslocamentos nos nós em y = 5 
   
        noscondess(end+1,:) = coords(i,1);

    end
end

% Distribuição das restrições pelos nós

for i = 1: size(noscondess, 1)

    no = noscondess(i);

    ixdof = (no-1)*2+1;    % Endereco do primeiro grau de liberdade deste no
    iydof = (no-1)*2+2;	    % Endereco do segundo grau de liberdade deste no

    Kr(iydof,iydof) = boom;     % Restringir graus de liberdade v (verticais)
    Kr(ixdof,ixdof) = boom;     % Restringir graus de liberdade u (horizontais)

    fr(iydof,1) = 0*boom;
    fr(ixdof,1) = 0*boom;

end

%------------------------------------------------------------------------------------------------------------------

% Arquivo da matriz global com condições naturais e essencias aplicadas

kr1 = Kr;


%------------------------------------------------------------------------------------------------------------------
                %% Resolução do sistema
%------------------------------------------------------------------------------------------------------------------

Kr = sparse(Kr);	    % Pode usar-se modo sparso para a matriz global

% Solução do sistema modificado por backslash

u = Kr\fr;

% Deslocamento máximo e deslocamento absoluto

umx = max(abs(u));

us = sqrt(u(1:2:Neqs,1).^2 + u(2:2:Neqs,1).^2);


% Verificacao do equilibrio estatico e cálculo das reações nos apoios	

R = Kg*u-fg;



%------------------------------------------------------------------------------------------------------------------
                %% Calcular Strain & Stresses
%------------------------------------------------------------------------------------------------------------------

% Vetores de tensões de von Mises e tensões normais

svms = [];
snorms = [];

%------------------------------------------------------------------------------------------------------------------

% Cálculo das tensões para cada elemento 

for i = 1:Nelt

    no1 = conects(i,2);
    no2 = conects(i,3);
    no3 = conects(i,4);

    Delta=zeros(6,1);


    % Endereços dos graus de liberdade dos 3 nós locais de cada triângulo 

    idof1 = (no1-1)*2+1;
    idof2 = (no2-1)*2+1;
    idof3 = (no3-1)*2+1;


    % Vetor com os graus de liberdade do nós do elemento

    edofs = [idof1 idof1+1 idof2 idof2+1 idof3 idof3+1];

    Delta = u(edofs);	% definir vector local de dofs ja calculados


    % Stresses e Strain calculado para cada elemento (A função foi retirada do
    % código Matlab do Professor Leonel Fernandes)
    
    [Stress Strain] = Stress_Strain_CST (x(no1),y(no1),x(no2),y(no2),x(no3),...
        y(no3),C,Delta);


    % Cálculo da tensão de von Mises para o caso de tensão plana

    I1= Stress(1)+Stress(2) ;               % Primeiro invariante
    I2= Stress(1)*Stress(2)-Stress(3)^2;  % Segundo invariante
    VMS = sqrt(I1*I1-3*I2);                  % von Mises Stress 

    % arquivo das tensões normais e de von Mises de todos os elementos num vetor

    svms(end+1, :) = VMS;
    snorms(end+1, :) = Stress(2);

end

%------------------------------------------------------------------------------------------------------------------

% Valores maximos e minimos das tensões

Tensao_Von_Mises_max = max(svms);
Tensao_Von_Mises_min = min(svms);
Tensao_Normal_Max = max(snorms);


%------------------------------------------------------------------------------------------------------------------
                %% Representação da malha pós-deformação
%------------------------------------------------------------------------------------------------------------------

% Coordenadas dos nós deformados

x1 = x + u(1:2:Neqs,1)*500;   % Factor de escala apropriado aqui
y1 = y + u(2:2:Neqs,1)*500;   % Factor de escala apropriado aqui

%------------------------------------------------------------------------------------------------------------------

% Representação gráfica

FIG10= figure (10);

subplot(1,2,1)

for i=1:Nelt

    no1=conects(i,2);
    no2=conects(i,3);
    no3=conects(i,4);
   
    edofs=[no1 no2 no3 no1];

    fill (x1(edofs),y1(edofs),us(edofs));hold on

    plot(x1(edofs),y1(edofs),'Color', [0.25, 0.25, 0.25]);hold on

end

plot(x1,y1,'.k');hold off

grid on

axis([-10, 10, -10, 10])

title('Solução: Deslocamentos')

%------------------------------------------------------------------------------------------------------------------

% Representaçao gráfica das tensões de von Mises

subplot(1,2,2)

for i=1:Nelt

    no1=conects(i,2);
    no2=conects(i,3);
    no3=conects(i,4);
   
    edofs=[no1 no2 no3 no1];
    tnsvms= [svms(i),svms(i), svms(i), svms(i) ];

    fill (x1(edofs),y1(edofs),tnsvms);hold on
    plot(x1(edofs),y1(edofs),'Color', [0.25, 0.25, 0.25]);hold on

end

grid on

axis([-10, 10, -10, 10])

title('Solução: Tensões de Von Mises')

%------------------------------------------------------------------------------------------------------------------

% Valores finais do exercício

com1 = [num2str(umx, '%e'),' mm'];
com2 = [num2str(Tensao_Von_Mises_max),' MPa'];

disp('Deslocamento Máximo'), disp(com1)
fprintf('\n')

disp('Tensão de Von Mises Máxima'), disp(com2) 
fprintf('\n')
