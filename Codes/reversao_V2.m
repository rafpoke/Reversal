%% Rotina de reversão
% A rotina de revesão cria um vetor de cargas simplificadas, as aplica na
% rotina de elementos finitos utilizada em estruturas
% A eficiência é definida como a sustentação adicional gerada pela
% superficie defletida deformada dividido pela sustentação adicional gerada
% pela superficie defletida rigida
% Para utilizar essa função, 
% *altere o load para o arquivo que você estiver usando como struct completo, 
% *altere os parametros do aileron de entrada
% *altere as funções criaGeometria, criaMateriais e criaLongarina para se
% adapatar para a sua realidade
%                               Boa sorte [Rafael, 2020]

clear all

load('o_escolhido.mat');  %% struct aviao completo a ser escolhido
ps = perfil_structure;


%% Parametros de entradas
Velocidades = 2:4:102;

Cl_d = 0.0637; % coeficiente de sustentação (derivada em relação à deflexão do aileron
fc_ail = 0.1258; % fator de corda do aileron 
fb_ail = 0.2745; % fator de envergadura do aileron (assume-se que o aileron existe de (1-fb)*b até a ponta de asa
pts_aileron = 5; % quantos pontos de pressão serão gerados no aileron (quantas conexões entre o aileron e a asa;mínimo 2)
S_ail = 0.021;   % area do aileron
da = 20;        % deflexão máxima do aileron


%% Cria o reversão struct

reversal = struct;
reversal.v = Velocidades;     %vetor das velocidades que serão testadas
reversal.efic = [];           %vetor com a eficiencia do aileron em cada velocidade
reversal.zero_state = [];     %vetor com a sustentação sem deflexão de profundor
reversal.rigid = [];          %vetor com a sustentação que a asa rígida geraria
reversal.flex = [];           %vetor com a sustentação gerada pela asa flexivel

%% a rotina

it_counter = 0;               % pra ver como está evoluindo a rotina 

for v=1:1:length(Velocidades)
    it_counter = it_counter +1
        % Cria cargas simplificadas 
        
        for i = 1:1:length(aviao.cl0_y1)
            alfa(i) = aviao.alfa_trim;

            F{1,1}(i,3) = (aviao.cl0_y1(i)+aviao.cl_alfa_y1(i)*deg2rad(alfa(i)))*const.rho*Velocidades(v)^2*aviao.b*aviao.c_y1(i)/(4*length(aviao.cl0_y1));
            F{1,1}(i,2) = 0;
            F{1,1}(i,1) = 0;

            P{1,1}(i,1) = 0.25*aviao.c_y1(i);
            P{1,1}(i,2) = 0.5*aviao.b*(i-1/2)/length(aviao.cl0_y1);
            P{1,1}(i,3) = 0.2585;
        end
        
        reversal.zero_state(v) = sum(F{1,1}(:,3));
        
        for j=1:1: pts_aileron
            F{1,1}(i+j,3) = (Cl_d*da)*const.rho*Velocidades(v)^2*S_ail/(2*pts_aileron);
            F{1,1}(i+j,2) = 0;
            F{1,1}(i+j,1) = 0;

            P{1,1}(i+j,1) = (1-fc_ail)*aviao.c_y1(i-pts_aileron+j);
            P{1,1}(i+j,2) = 0.5*aviao.b*(i-pts_aileron+j)/length(aviao.cl0_y1);
            P{1,1}(i+j,3) = 0.2585; %pra asa de cima adiconar a distancia entre asas
        end

        cargas_simples.F = F;
        cargas_simples.P = P;

        
        reversal.rigid(v) = sum(F{1,1}(:,3));
        
        % Parte de Estruturas

        % ============================= Struct geom ============================= %
        [ geom_Asa ] = ETT_criaGeom2( aviao , ps );      % structs com informações geométricas das supercícies aerodinâmicas
        [ geom_Asa ] = ETT_criaEH(aviao, ps, geom_Asa);
        % =========================== Struct materiais ========================== %
        [ materiais ] = ETT_criaMateriais2020();                                        % struct com propriedades dos materiais
        % =========================== Struct Longarina ========================== %
        % [ L_w ] = ETT_criaLongarinaRelatorio( materiais );                                   % struct com dados das longarinas
        load('Long_RafaEscolhidaHW.mat');


        % elementos finitos

        L_w = ETT_PontosLongarinaGeral2( geom_Asa , L_w, ps );

        L_w = ETT_PontosTirante( geom_Asa , L_w );

        L_w = ETT_PropriedadesLongarina( L_w );

        L_w = ETT_Main_FEM( L_w , cargas_simples );

        % repetir
        for i = 1:1:length(alfa)
            beta(i) = alfa(i) + rad2deg(L_w(1).u(5, int16(length(L_w(1).u)*i/length(alfa))));
        end
        
        % calcula cargas atualizadas
        for i = 1:1:length(aviao.cl0_y1)
           carga_z_updated(i) = (aviao.cl0_y1(i)+aviao.cl_alfa_y1(i)*deg2rad(beta(i)))*const.rho*Velocidades(v)^2*aviao.b*aviao.c_y1(i)/(4*length(aviao.cl0_y1));
        end    
        
        reversal.flex(v) = sum(carga_z_updated)+((Cl_d*da)*const.rho*Velocidades(v)^2*S_ail/2);
        reversal.efic(v) = (reversal.flex(v) - reversal.zero_state(v))/(reversal.rigid(v)-reversal.zero_state(v));
        
        clear geom_Asa
        clear materiais
        clear L_w

        
        if reversal.efic(v)<=0
            reversal.v = Velocidades(1:v);
            break
        end
end

%% plotar o final

figure(1)
hold on
plot(reversal.v, reversal.efic,'b')
plot([22,22],[-100,100]);
hold off

figure(2)
hold on
plot(reversal.v, reversal.flex - reversal.zero_state, '*' )
plot(reversal.v, reversal.rigid - reversal.zero_state, '--')
hold off