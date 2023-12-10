
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
% - Simulação do caso simples (quadrado sem entalhe 10mmx10mm) através de MEF por Triângulo de 6 Nós
% - Partes do código foram retiradas do código proporconado nas aulas do
% professor Leonel Fernandes
%------------------------------------------------------------------------------------------------------------------


%------------------------------------------------------------------------------------------------------------------
                %% Comandos de Limpeza
%------------------------------------------------------------------------------------------------------------------

clear all
close all

%------------------------------------------------------------------------------------------------------------------
                %% Especificações do Elemento
%------------------------------------------------------------------------------------------------------------------

%  Dimensões

h = 10;
esp = 1;

% Propriedades do material

E = 70*10^3;
niu = 0.3;
G = E/(2*(1+niu));

% Relação constitutiva

cons=E/(1 - niu*niu);
C=cons*[1 niu 0;niu 1 0;0 0 (1-niu)/2];



%------------------------------------------------------------------------------------------------------------------
                %% Importação da malha do NX através do ficheiro txt
%------------------------------------------------------------------------------------------------------------------

% Importante referir que a origem do referencial no modelo da peça encontra
% se no canto inferior esquerdo da mesma para todos os casos. Este programa
% automaticamente corrige coordenadas que não obedeçam o esquema


% Importação do ficheiro txt



% Caso a orientação dos nós esteja trocada alterar o valor lógico da
% variavel para TRUE

inver = false;


% Para o cálculo das caracteristicas especificas e resultados em pontos
% objetivos de estudo foram selecionados 3 pontos, indicados no relatório,
% pontos 1, 2 e 3. No caso da geometria simples o ponto 2 não é considerado
% Caso a geometria proposta para correr neste código não apresente o pinto
% 2 alterar o valor lógico da variavel seguinte para FALSE

ponto_2_on = false;

fileID0 = fopen("NOSfem1T6GeomSimples.txt", "r");
fileID1 = fopen("ELEMfem1T6GeomSimples.txt", "r");
formatSpec = '%c'; '%d';


% Leitura dos valores das coordenadas dos nós


nos = splitlines(fscanf(fileID0, formatSpec));
aux = [];

for i = 11:6:size(nos, 1)

    aux = [aux, nos(i)];

end

aux=aux';
colunas = split(aux);

coord1 = [];
coord2 = [];
coordx = [];
coordy = [];
coordout = [];

% O sistema de coordenadas é dado em coordout

for i=1:1:size(colunas, 1)
    coord1 = [coord1; colunas(i, 5)];
    coordx = [coordx; str2double(coord1(i))];
    coord2 = [coord2; colunas(i, 6)];
    coordy = [coordy; str2double(coord2(i))];
    coordout = [coordout; i, coordx(i), coordy(i)];
end
  

% Leitura dos valores dos nós correspondentes a cada elemento


elementos = splitlines(fscanf(fileID1, formatSpec));

aux2 = [];

for i = 17:1:(size(elementos, 1) - 4)

    aux2 = [aux2, elementos(i)];

end

aux2=aux2';

coluna2 = [];
coluna2 = split(aux2);

connod1 = [];
connod2 = [];
connod3 = [];
connod4 = [];
connod5 = [];
connod6 = [];
connodout = [];

% O conjunto dos nós referente a cada elemento é dado em connodout

for i = 1:1:size(coluna2, 1)

    connod1 = [connod1; str2double(coluna2(i, 12))];
    connod2 = [connod2; str2double(coluna2(i, 13))];
    connod3 = [connod3; str2double(coluna2(i, 14))];
    connod4 = [connod4; str2double(coluna2(i, 15))];
    connod5 = [connod5; str2double(coluna2(i, 16))];
    connod6 = [connod6; str2double(coluna2(i, 17))];
    connodout = [connodout; i, connod1(i), connod2(i), connod3(i), connod4(i), connod5(i), connod6(i)];

end



%------------------------------------------------------------------------------------------------------------------
                %% Reorganização para acertar a orientação dos nós e correção do referencial
%------------------------------------------------------------------------------------------------------------------


% De seguida correm linhas de comando para reorganizar a matriz de
% conectividades de modo a arranjar erros de sentido dos nós. (o sentido considerado pelo
% código é o anti-horário

if inver == true

    connodout(:, [2,4]) = fliplr(connodout(:, [2,4]));

    connodout(:, [5,6]) = fliplr(connodout(:, [5,6]));

end

% O código seguinte orienta a figura para um referencial de origem no nó do
% canto inferior esquerdo da componente

size(coordout,1); 

coordmaxx = 0;
coordmaxy = 0;


for i = 1:size(coordout,1)
    
    coordxteste = abs(coordout(i, 2));
    coordyteste = abs(coordout(i, 3));

    if coordxteste > coordmaxx

        coordmaxx = coordxteste;
    end

    if coordyteste > coordmaxy

        coordmaxy = coordyteste;
    end

end

move_dir = 10 - coordmaxx;
move_up = 10 - coordmaxy;

for i = 1:size(coordout,1)
    
   coordout(i, 2) = coordout(i, 2) + move_dir;
   coordout(i, 3) = coordout(i, 3) + move_up;

end



% Designação dos valores importados para variaveis adaptadas ao código

x= [0; 10; 10; 0];
y= [0; 0; 10; 10];

tri6 = connodout(:, [2 3 4 5 6 7]);
coords = coordout;

x = coords(:,2);
y = coords(:,3);


%------------------------------------------------------------------------------------------------------------------
                %% Representação gráfica do elemento 
%------------------------------------------------------------------------------------------------------------------

FIG3 = figure(3)
Nelt = size(tri6,1); 

%------------------------------------------------------------------------------------------------------------------

for i=1:Nelt

    no1 = tri6(i,1); 
    no2 = tri6(i,2); 
    no3 = tri6(i,3); 
    no4 = tri6(i,4); 
    no5 = tri6(i,5); 
    no6 = tri6(i,6); 

    edofs = [tri6(i,1) tri6(i,4) tri6(i,2) tri6(i,5) tri6(i,3) tri6(i,6)]; 

    fill (x(edofs),y(edofs),'green');hold on
    plot(x(edofs),y(edofs),'black');hold on

    % Possibilidade de representar nos centros o número com o elemento

     % cenx= (x(no1) + x(no2) + x(no3))/3;
     % ceny= (y(no1) + y(no2) + y(no3))/3;
     % text(cenx, ceny, num2str(i)); hold on

end

% plot da geometria

plot(x,y,'.'); hold on
axis([-5, 15, -5, 15]);
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
% função Elem_ElasTRI6 e assemblagem da matriz global

for i = 1:Nelt
    
    no1 = tri6(i,1); 
    no2 = tri6(i,2); 
    no3 = tri6(i,3); 
    no4 = tri6(i,4); 
    no5 = tri6(i,5); 
    no6 = tri6(i,6); 

    
    edofs = [no1 no2 no3 no4 no5 no6];
    XN(1:6,1) = x(edofs);
    XN(1:6,2) = y(edofs);
    
    % Endereços na matriz global dos graus de liberdade dos 6 nós locais de
    % cada triângulo 
    
    idof1 = (no1-1)*2+1;
    idof2 = (no2-1)*2+1;
    idof3 = (no3-1)*2+1;
    idof4 = (no4-1)*2+1;
    idof5 = (no5-1)*2+1;
    idof6 = (no6-1)*2+1;
    
    % Vetor com os graus de liberdade do nós do elemento

    edofs = [idof1 idof1+1 idof2 idof2+1 idof3 idof3+1 idof4 idof4+1 idof5 ...
         idof5+1 idof6 idof6+1]; 

    % Matriz local do elemento e vetor força (A função foi retirada do
    % código Matlab do Professor Leonel Fernandes)

    [Ke fe] = Elem_ElasTRI6 (XN,C,esp,0);

    % Assemblagem do elemento associando a matriz local na matriz global e
    % o vetor forças local no global

    Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;
    fg(edofs,1)= fg(edofs,1) + fe;

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
 

noscondat= [];

for i= 1:Nnds

    if y(i) == 10

        noscondat(end+1,:) = coords(i,1);

    end

end

% Número de nós com condições naturais e matriz das conectividades

Nnoscondat = size(noscondat, 1);
connodout2 = connodout(:, [2 3 4 5 6 7]);

% Cálculos e códigos para as distribuições da tração nas faces quadráticas
% dos triângulos

comp1 = h/(Nnds - 1); % comprimento entre cada nó

noscondatinter = [];    % Nós com condições naturais em pontos intermédios dos nós principais
noscondatprin = [];     %Nós com condiçoes naturais em cantos do triângulo (nós principais)

% Código para assimilar estas duas listas dentro dos nós com condições
% naturais para masi tarde aplicar as trações segundo a fórmula do integral
% numérco para o elemento unidimensional quadrático


for i = 1:Nnoscondat

    no = noscondat(i, 1);
    
    [row, col] = find(connodout2 == no);

    if col < 4

        noscondatprin(end+1,:) = no;

    else 

        noscondatinter(end+1,:) = no;

    end

end

% Aplicação da carga distribuida

qx = 0;
qy = 2;
fpnx = (qx*h/(3*Nnoscondat - 3));
fpny = (qy*h/(3*Nnoscondat - 3));

for i = 1:size(noscondatprin)

    nop = noscondatprin(i);

    if coords(nop, 2) == 0

        % disp('Canto esquerdo encontrado')

        idof = (nop-1)*2 + 2;
        fg(idof,1) = fg(idof,1) + fpny;

    
    elseif coords(nop, 2) == 10

        % disp('Canto direito encontrado')
        
        idof = (nop-1)*2 + 2;

        fg(idof,1) = fg(idof,1) + fpny;
        
    else

        iydof = (nop-1)*2 + 2;       % Endereço do segundo grau de liberdade deste nó
        fg(iydof,1) = fg(iydof,1) + fpny*2;     % Força aplicada no nó

    end

end



for i = 1:size(noscondatinter)

    noi = noscondatinter(i);

    if coords(noi, 2) == 0

        disp('Detetado nó intermedio num canto! Erro')
        
    elseif coords(noi, 2) == 10

        disp('Detetado nó intermedio num canto! Erro')
      
    else

        iydof = (noi-1)*2 + 2;       % Endereço do segundo grau de liberdade deste nó
        fg(iydof,1) = fg(iydof,1) + fpny*4;        % Força aplicada no nó

    end

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

for i= 1:Nnds

    if y(i) == 0        % Aplicação da restrição dos deslocamentos nos nós em y = 5 

        noscondess(end+1,:) = coords(i,1);

    end

end

% Distribuição das restrições pelos nós

for i= 1: size(noscondess, 1)

    no= noscondess(i);

    ixdof = (no-1)*2+1;    % Endereco do primeiro grau de liberdade deste no
    iydof = (no-1)*2+2;    % Endereco do segundo grau de liberdade deste no	

    Kr(iydof,iydof) = boom;     % Restringir graus de liberdade v (verticais)
    Kr(ixdof,ixdof) = boom;     % Restringir graus de liberdade u (horizontais)
    
    fr(iydof,1)= 0*boom;
    fr(ixdof,1)= 0*boom;

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

% Energia de deformação elástica


Energ_Def_El = (1/2) * u' * fg;


%------------------------------------------------------------------------------------------------------------------
                %% Calcular Strain & Stresses
%------------------------------------------------------------------------------------------------------------------

% Vetores de tensões de von Mises e tensões normais

svms = [];
snorms = [];
Stress_Total = [];

%------------------------------------------------------------------------------------------------------------------

% Cálculo das tensões para cada elemento 

for i = 1:Nelt

    Delta = zeros(12,1);

     % Endereços dos graus de liberdade dos 6 nós locais de cada triângulo 

    no1 = tri6(i,1); 
    no2 = tri6(i,2); 
    no3 = tri6(i,3); 
    no4 = tri6(i,4); 
    no5 = tri6(i,5); 
    no6 = tri6(i,6); 

    % Vetor com os graus de liberdade do nós do elemento

    edofs =[no1 no2 no3 no4 no5 no6];


    % Coordenadas do elemento

    XN(1:6,1)=x(edofs);
    XN(1:6,2)=y(edofs);

     % Endereços dos graus de liberdade dos 3 nós locais de cada triângulo 

    idof1 = (no1-1)*2+1;
    idof2 = (no2-1)*2+1;
    idof3 = (no3-1)*2+1;
    idof4 = (no4-1)*2+1;
    idof5 = (no5-1)*2+1;
    idof6 = (no6-1)*2+1;

    % Vetor com os graus de liberdade do nós do elemento

    edofs = [idof1 idof1+1 idof2 idof2+1 idof3 idof3+1 idof4 idof4+1 idof5 ...
         idof5+1 idof6 idof6+1];


    Delta = u(edofs);	% definir vector local de dofs ja determinados

    % Centroide do elemento

    csi = 1/3;
    eta = 1/3;

     % Stresses e Strain calculado para cada elemento (A função foi retirada do
     % código Matlab do Professor Leonel Fernandes)

     [Stress Strain] = Stress_Strain_TRI6 (XN,C,Delta,csi,eta);


    % Cálculo da tensão de von Mises para o caso de tensão plana

    I1= Stress(1)+Stress(2) ;               % Primeiro invariante
    I2= Stress(1)*Stress(2)-Stress(3)^2;  % Segundo invariante
    VMS = sqrt(I1*I1-3*I2);                  % von Mises Stress 

     % arquivo da tensão de von Mises de todos os elementos num vetor
    
     svms(end+1, :) = VMS;
     snorms(end+1, :) = Stress(2);
     Stress_Total = [Stress_Total; Stress'];

end 

%------------------------------------------------------------------------------------------------------------------

% Valores maximos e minimos das tensões

Tensao_Von_Mises_max = max(svms);
Tensao_Von_Mises_min = min(svms);


%------------------------------------------------------------------------------------------------------------------

% Cálculos de pontos de interesse


% Indexes dos nós pretendidos


index1 = find(coords(:, 2) == 10 & coords(:, 3) == 10);
index3 = find(coords(:, 2) == 10 & coords(:, 3) == 0);

if ponto_2_on == true

    index2 = find(coords(:, 2) == 7.9852 & coords(:, 3) == 4.9884);
else 
    index2 = 0;
end

% elementos em que partilham os nós pretendidos

elem_no1 = find(connodout(:, 2) == index1 | connodout(:, 3) == index1 | connodout(:, 4) == index1);
if ponto_2_on == true
    elem_no2 = find(connodout(:, 2) == index2 | connodout(:, 3) == index2 | connodout(:, 4) == index2);
end
elem_no3 = find(connodout(:, 2) == index3 | connodout(:, 3) == index3 | connodout(:, 4) == index3);

%------------------------------------------------------------------------------------------------------------------

% Cálculo das medias das tensões de von Mises dos elementos á volta dos nós
% pretendidos

% Nó 1

vms_no1 = 0;

for i = 1:numel(elem_no1)

    vms_no1 = vms_no1 + svms(elem_no1(i));

end

vms_no1 = vms_no1/numel(elem_no1);

% Nó 2

if ponto_2_on == true

    vms_no2 = 0;

    for i = 1:numel(elem_no2)

        vms_no2 = vms_no2 + svms(elem_no2(i));

    end

vms_no2 = vms_no2/numel(elem_no2);

end

% Nó 3

vms_no3 = 0;

for i = 1:numel(elem_no3)

    vms_no3 = vms_no3 + svms(elem_no3(i));

end

vms_no3 = vms_no3/numel(elem_no3);

%------------------------------------------------------------------------------------------------------------------

% Deslocamentos em y de cada ponto

uy_no1 = u(2*(index1), 1);
if ponto_2_on == true
    uy_no2 = u(2*(index2), 1);
end
uy_no3 = u(2*(index3), 1);

%------------------------------------------------------------------------------------------------------------------

% Cálculo das medias das tensões de von Mises dos elementos á volta dos nós
% pretendidos

% Nó 1

sigmayy_no1 = 0;

for i = 1:numel(elem_no1)

    sigmayy_no1 = sigmayy_no1 + snorms(elem_no1(i));

end

sigmayy_no1 = abs(sigmayy_no1/numel(elem_no1));

% Nó 2

list_sigma_el = [];

if ponto_2_on == true

    sigmayy_no2 = 0;

    
    for i = 1:numel(elem_no2)

        sigmayy_no2 = sigmayy_no2 + snorms(elem_no2(i));
        list_sigma_el = [list_sigma_el, Stress_Total(elem_no2(i), :)'];           % Lista dos sigmas dos elementos em volta do nó 2 para o cálculo da energia de deformação elástica neste ponto

    end



    sigmayy_no2 = abs(sigmayy_no2/numel(elem_no2));

end

% Nó 3

sigmayy_no3 = 0;

for i = 1:numel(elem_no3)

    sigmayy_no3 = sigmayy_no3 + snorms(elem_no3(i));

end

sigmayy_no3 = abs(sigmayy_no3/numel(elem_no3));

%------------------------------------------------------------------------------------------------------------------

% Cálculo da energia de deformação elástica por unidade de volume do nó 2 através da média das
% energias de deformação dos elementos adjacentes a este nó. A fórmula foi
% retirada do livro "Theory of Elasticity", S. Timoshenko and J. N.
% Goodier. Pág 84 - (86)

% Cálculo da energia para cada elemento
if ponto_2_on == true
    Soma_Vo = 0;
    list_sigma_el = list_sigma_el';

    for i = 1:numel(elem_no2)

        Vo = (0.5)*(1/E)*( ((list_sigma_el(i, 1))^2) + ((list_sigma_el(i, 2))^2) ) - (niu/E)*(list_sigma_el(i, 1))*(list_sigma_el(i, 2)) + 0.5*(1/G)*(list_sigma_el(i, 3));
        Soma_Vo = Soma_Vo + Vo;

    end

    Energ_Def_El_no2 = Soma_Vo/numel(elem_no2);

end

%------------------------------------------------------------------------------------------------------------------
                %% Representação da malha pós-deformação
%------------------------------------------------------------------------------------------------------------------

% Coordenadas dos nós deformados

x1 = x + u(1:2:Neqs,1)*500;   % Factor de escala apropriado aqui
y1 = y + u(2:2:Neqs,1)*500;   % Factor de escala apropriado aqui

%------------------------------------------------------------------------------------------------------------------

% Representação gráfica

FIG11= figure (11);

subplot(1,2,1)

for i = 1:Nelt

    no1 = tri6(i,1); 
    no2 = tri6(i,2); 
    no3 = tri6(i,3); 
    no4 = tri6(i,4); 
    no5 = tri6(i,5); 
    no6 = tri6(i,6); 

    edofs=[tri6(i,1) tri6(i,4) tri6(i,2) tri6(i,5) tri6(i,3) tri6(i,6)];

    plot(x1(edofs),y1(edofs),'Color', [0.25, 0.25, 0.25]);hold on

    fill(x1(edofs), y1(edofs), us(edofs)); hold on

end


plot(x1,y1,'.k');hold off

axis([-5, 15, -5, 15])

grid on

title('Solução: Deslocamentos')

%------------------------------------------------------------------------------------------------------------------

% Representaçao gráfica das tensões de von Mises

subplot(1,2,2)

for i = 1:Nelt
    no1 = tri6(i,1); 
    no2 = tri6(i,2); 
    no3 = tri6(i,3); 
    no4 = tri6(i,4); 
    no5 = tri6(i,5); 
    no6 = tri6(i,6); 

    edofs=[tri6(i,1) tri6(i,4) tri6(i,2) tri6(i,5) tri6(i,3) tri6(i,6)];

    tnsvms= [svms(i), svms(i), svms(i), svms(i), svms(i), svms(i),];

    fill(x1(edofs), y1(edofs), tnsvms); hold on

    plot(x1(edofs),y1(edofs),'Color', [0.25, 0.25, 0.25]);hold on

end

axis([-5, 15, -5, 15])

grid on

title('Solução: Tensões de Von Mises')

%------------------------------------------------------------------------------------------------------------------

% Valores finais da simulação

com1 = [num2str(umx, '%e'),' mm'];
com2 = [num2str(Tensao_Von_Mises_max),' MPa'];
com3 = [num2str(Energ_Def_El, '%e'),' mJ'];

com1no1 = [num2str(vms_no1),' MPa'];
if ponto_2_on == true
    com1no2 = [num2str(vms_no2),' MPa'];
end
com1no3 = [num2str(vms_no3),' MPa'];

com2no1 = [num2str(uy_no1, '%e'),' mm'];
if ponto_2_on == true
    com2no2 = [num2str(uy_no2, '%e'),' mm'];
end
com2no3 = [num2str(uy_no3, '%e'),' mm'];

com3no1 = [num2str(sigmayy_no1),' MPa'];
if ponto_2_on == true
    com3no2 = [num2str(sigmayy_no2),' MPa'];
end
com3no3 = [num2str(sigmayy_no3),' MPa'];

if ponto_2_on == true
    com4no2 = [num2str(Energ_Def_El_no2, '%e'),' mJ/mm^3'];
end

fprintf('\n')
fprintf('\n')
disp('Deslocamento Máximo'), disp(com1)
fprintf('\n')

disp('Tensão de Von Mises Máxima'), disp(com2) 
fprintf('\n')

disp('Energia de Deformação Elástica'), disp(com3)
fprintf('\n')

fprintf('\n')
fprintf('\n')
disp('Valores para pontos especificos')
fprintf('\n')
disp('Para o nó 1')
fprintf('\n')
disp('Tensão de Von Mises'), disp(com1no1)
fprintf('\n')
disp('Deslocamento em y'), disp(com2no1)
fprintf('\n')
disp('Tensão em y'), disp(com3no1)
fprintf('\n')

if ponto_2_on == true

    fprintf('\n')
    fprintf('\n')
    disp('Para o nó 2')
    fprintf('\n')
    disp('Tensão de Von Mises'), disp(com1no2)
    fprintf('\n')
    disp('Deslocamento em y'), disp(com2no2)
    fprintf('\n')
    disp('Tensão em y'), disp(com3no2)
    fprintf('\n')
    disp('Energia de Deformação Elástica'), disp(com4no2)
    fprintf('\n')
    
end

fprintf('\n')
fprintf('\n')
if ponto_2_on == true
    disp('Para o nó 3')
else
    disp('Para o nó 2')
end
fprintf('\n')
disp('Tensão de Von Mises'), disp(com1no3)
fprintf('\n')
disp('Deslocamento em y'), disp(com2no3)
fprintf('\n')
disp('Tensão em y'), disp(com3no3)
fprintf('\n')