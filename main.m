
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
% - Menu dos vários programas
% - Partes do código dos programas foram retiradas do código proporconado nas aulas do
% professor Leonel Fernandes
%------------------------------------------------------------------------------------------------------------------


%------------------------------------------------------------------------------------------------------------------
                %% Menu
%------------------------------------------------------------------------------------------------------------------

disp('Dentro das simulações  ')
a = input(['Qual o programa a correr? \n1: Simulação da malha com triangulos T3 no caso da geometria simples  \n2: Simulação da malha proporconada pelo NX com triangulos T3 \n3: Simulação da malha proporconada pelo NX com triangulos T6 \n4: Sair \n \n'], 's');

if a == '1'

    % Simulação da malha com triangulos T3 bo Caso da geometria simples

    SIM_T3_CasoSimples

    drawnow

    input('Enter para continuar')

    close all

    main

elseif a== '2'

    % Simulação da malha proporconada pelo NX com triangulos T3

    SIM_T3_NX

    drawnow

    input('Enter para continuar')

    close all

    main
elseif a== '3'

    % Simulação da malha com triangulos T6

    SIM_T6_NX

    drawnow

    input('Enter para continuar')

    close all

    main

elseif a == '4'

    % Sair do menu

    close all

    disp('Saiu do menu')

else

    disp('Letra errada!')
    main

end



