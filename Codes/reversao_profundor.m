%% Rotina de revers�o de empenagem
% A rotina de reves�o cria um vetor de cargas simplificadas, as aplica na
% rotina de elementos finitos utilizada em estruturas
% A efici�ncia � definida como a sustenta��o adicional gerada pela
% superficie defletida deformada dividido pela sustenta��o adicional gerada
% pela superficie defletida rigida
% Para utilizar essa fun��o, 
% *altere o load para o arquivo que voc� estiver usando como struct completo, 
% *altere os parametros do aileron de entrada
% *altere as fun��es criaGeometria, criaMateriais e criaLongarina para se
% adapatar para a sua realidade
% Tive problemas pra repetir a fun��o pra empenagem horizontal
% criei uma igual mas diferente
%                               Boa sorte [Rafael, 2020]

clear all
close all

load('o_escolhido.mat');  %% struct aviao completo a ser escolhido
ps = perfil_structure;


%% Parametros de entradas
Velocidades = 2:4:102;

Cl_de = -0.0115; % coeficiente de sustenta��o (derivada em rela��o � deflex�o do aileron
c_prof = 0.1805;
fc_prof = 0.2963*c_prof; % fator de corda do profundor
b_prof = 1.9699/2;
fb_prof = 0.7928*b_prof; % fator de envergadura do profundor (assume-se que o aileron existe de (1-fb)*b at� a ponta de asa
pts_profundor = 10; % quantos pontos de press�o ser�o gerados no aileron (quantas conex�es entre o aileron e a asa;m�nimo 2)
S_prof = 0.3387;   % area do aileron
da = 20;        % deflex�o m�xima do profundor

rho = const.rho;


%% Cria o revers�o struct

reversal = struct;
reversal.v = Velocidades;     %vetor das velocidades que ser�o testadas
reversal.efic = [];           %vetor com a eficiencia do aileron em cada velocidade
reversal.zero_state = [];     %vetor com a sustenta��o sem deflex�o de profundor
reversal.rigid = [];          %vetor com a sustenta��o que a asa r�gida geraria
reversal.flex = [];           %vetor com a sustenta��o gerada pela asa flexivel

%% a rotina

it_counter = 0;               % pra ver como est� evoluindo a rotina :D

for v=1:1:length(Velocidades)
    it_counter = it_counter +1
        % Cria cargas simplificadas 
                
        for j=1:1: pts_profundor
            F{1,1}(j,3) = (-Cl_de*da)*const.rho*Velocidades(v)^2*S_prof/(2*pts_profundor);
            F{1,1}(j,2) = 0;
            F{1,1}(j,1) = 0;

            P{1,1}(j,1) = 0.46+(c_prof - fc_prof);
            P{1,1}(j,2) = (b_prof - fb_prof*(pts_profundor-j)/pts_profundor);
            P{1,1}(j,3) = 0.65; %pra asa de cima adiconar a distancia entre asas
        end

        cargas_simples.F = F;
        cargas_simples.P = P;

        
        reversal.rigid(v) = sum(F{1,1}(:,3));
        
        % Parte de Estruturas

        % ============================= Struct geom ============================= %
        [ geom_Asa ] = ETT_criaGeom2( aviao , ps );      % structs com informa��es geom�tricas das superc�cies aerodin�micas
        [ geom_Asa ] = ETT_criaEH(aviao, ps, geom_Asa);

        % =========================== Struct materiais ========================== %
        [ materiais ] = ETT_criaMateriais2020();                                        % struct com propriedades dos materiais
        % =========================== Struct Longarina ========================== %
        % [ L_w ] = ETT_criaLongarinaRelatorio( materiais );                                   % struct com dados das longarinas
        load('L_w_EH.mat');



        % elementos finitos

        L_w = ETT_PontosLongarinaGeral2( geom_Asa , L_w, ps );

        L_w = ETT_PontosTirante( geom_Asa , L_w );

        L_w = ETT_PropriedadesLongarina( L_w );

        L_w = ETT_Main_FEM( L_w , cargas_simples );

        % repetir
        for i = 1:1:10
            beta(i) = rad2deg(L_w(5).u(5, int16(length(L_w(1).u)*i/10)));
        end
        
        % calcula cargas atualizadas
        for i = 1:1:10
           carga_z_updated(i) = 2*pi*deg2rad(beta(i))*const.rho*Velocidades(v)^2*S_prof/(2*10);
        end    
        
        reversal.flex(v) = sum(carga_z_updated)+((-Cl_de*da)*const.rho*Velocidades(v)^2*S_prof/2);
        reversal.efic(v) = (reversal.flex(v))/(reversal.rigid(v));
        
        clear geom_Asa
        clear materiais


        
        if reversal.efic(v)<=0
            reversal.v = Velocidades(1:v);
            break
        end
end

reversal.zero_state = 0;

%% plotar o final

figure(1)
hold on
plot(reversal.v, reversal.efic,'b')
plot([22,22],[-1,2])
hold off

figure(2)
hold on
plot(reversal.v, reversal.flex, '*' )
plot(reversal.v, reversal.rigid, '--')
hold off