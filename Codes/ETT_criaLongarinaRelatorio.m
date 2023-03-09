%% Fun��o ETT_criaLongarina %%
% ============================== DESCRI��O ============================== %
% Essa fun��o parametriza completamente a geometria de longarinas com
% m�ltiplas se��es. Essa parametriza��o � �til no caso de integra��o com
% uma rotina de otimiza��o. Os par�metros podem ser de dois tipos: globais
% ou espec�ficos. Par�metros globais s�o utilizados para descrever qualquer
% tipo de longarina (independentemente da geometria da se��o). Par�metros
% espec�ficos descrevem uma geometria espec�fica (como circular, retangular
% etc). Para auxiliar o entendimento da fun��o � recomendado alterar
% par�metros e utilizar a fun��o plotLongarina() para melhor visualiza��o.

% =============================== OUTPUTS =============================== %
% - Longarina - struct [1xn] onde n � o n�mero de longarinas parametrizadas
% da aeronave.

% =========================== PARAMETRIZA��O ============================ %

% =============================== GLOBAL ================================ %
% - config - cell array [1xk] (onde k � o n�mero de se��es da longarina) de
% strings. As string simbolizam o formato da se��o da longarina. Formatos
% poss�veis: 'circular' (programado), 'retangular' (n�o-programado), 'I'
% (n�o-programado), 'C' (n�o-programado), 'BA' (n�o-programado), 'tubo'
% (n�o-programado).
% Ex: Longarina(1).config = {'circular','retangular'}

% - material - cell array [1xk] (onde k � o n�mero de se��es da longarina)
% de structs. As structs indicam o material do qual � feita a se��o de
% longarina. Materiais poss�veis: carbono (programado), fibra_vidro
% (n�o-programado), dyneema (n�o-programado), kevlar_carbono 
% (n�o-programado), balsa (n�o-programado), plywood (n�o-programado).
% Ex: Longarina(3).material = {materiais.carbono,materiais.balsa}

% - cc - cell array [nx2] (onde n � o n�mero de condi��es de contorno da
% longarina). Cada linha indica o n�mero do n� em que se aplica a condi��o
% de contorno e os graus de liberdade restringidos de 1 a 6 (u v w tx ty
% tz).

% - fb - vetor [1x(k-1)] (onde k � o n�mero de se��es da longarina) da
% fra��o da envergadura restante ocupada pela se��o da longarina. Por
% exemplo em uma asa de semi_envergadura b=1 e fb = [0.5 0.5] ter�amos uma
% longarina com tr�s se��es de envergaduras [0.5 0.25 0.25]; no caso de fb
% = [0.4 0.3] ter�amos envergaduras [0.4 0.18 0.42].
% Ex: Longarina(1).fb = [0.4] ou Longarina(1).fb = [] (no caso de uma se��o
% �nica de longarina)

% - fc - vetor [1x(k+1)] (onde k � o n�mero de se��es da longarina) que 
% indicam a fra��o de corda em que est�o posicionados os pontos de
% troca de se��o da longarina.
% Ex: Longarina(2).fc = [0.5 0.4 0.25] (nesse caso a se��o 1 teria
% extremidades em 50% de corda e 40% de corda, enquanto a se��o 2 teria
% extremidades em 40% de corda e 25% de corda)

% - fh - vetor [1x(k+1)] (onde  k � o n�mero de se��es da longarina) que
% indica o qu�o afastado em z do centro do perfil est� o centro da se��o
% da longarina. Esse valor � uma fra��o da dist�ncia em z da parte superior
% da asa � parte inferior que pode variar entre [-0.5,+0.5].
% Ex: Longarina(1).fh = [0 -0.05 0.1]

% ============================== CIRCULAR =============================== %
% - fr - vetor [1xk] (onde  k � o n�mero de se��es da longarina) que indica
% a fra��o do raio m�ximo da se��o que ser� utlizado. Por exemplo, se a asa
% comportaria uma se��o j de longarina circular com raio 15mm e fr(j) = 0.5
% ent�o a se��o j da longarina ter� raio de 7.5mm.
% Ex: Longarina(2).fr = [0.7 0.6]

% - f_voltas - vetor [1xk] (onde  k � o n�mero de se��es da longarina) que
% indica a fra��o de envergadura da se��o em que ocorrer� uma mudan�a de
% n�mero de voltas de laminado. Por exemplo, se a envergadura da se��o j �
% b=0.6 e f_voltas(j)=0.4 ent�o ocorrer� uma troca de n�mero de voltas em
% y=0.24.
% Ex: Longarina(2).f_voltas = [0.2 nan 0.8] (no caso das se��es [1,3] serem
% circulares e a se��o [2] n�o)

% - n_voltas - matriz [kx2] (onde  k � o n�mero de se��es da longarina) em
% que cada linha indica o n�mero de voltas de laminado de cada se��o de
% longarina.
% Ex: Longarina(3).n_voltas = [3 2;nan nan;2 2] (nesse caso na primeira
% se��o a longarina inicialmente tem 3 voltas e depois tem 2. A segunda
% se��o n�o � circular e a terceira tem 2 voltas ao longo de toda a
% envergadura)

% ============================= RETANGULAR ============================== %
% phi [deg] - vetor [1x(k+1)] (onde k � o n�mero de se��es da longarina)
% que indica a rela��o entre a alma e a base (phi=atand(alma/base)) em cada
% troca de se��o da longarina.
% Ex: Longarina(1) = [45 45 45] (nesse caso a longarina tem 2 se��es, todas
% com alma igual � base).

% e_alma [m] - vetor [1xk] (onde k � o n�mero de se��es da longarina) que
% indica a espessura da alma da longarina retangular em cada se��o. A alma 
% � a lateral da longarina.

% e_base [m] - vetor [1xk] (onde k � o n�mero de se��es da longarina) que
% indica a espessura da base da longarina retangular em cada se��o. A base 
% � a parte superior e inferior da longarina.

% f_dim - vetor [1x(k+1)] (onde k � o n�mero de se��es da longarina) que 
% indica a fra��o da dimens�o m�xima dispon�vel em uma troca de se��o que
% a longarina ocupar�.
% Ex: Longarina(1).f_dim = [1 0.9] (nesse caso a longarina possui apenas
% uma se��o. Se na ra�z a dimens�o m�xima dispon�vel para a alma fosse 40mm
% e na ponta fosse 10mm, as dimens�es reais da alma seriam 40mm na ra�z e
% 9mm na ponta.

% ============================== SANDUICHE ============================== %
% D_I [N.m] - vetor [1x2] que indica a rigidez � flex�o por unidade de
% largura (EI/b) da se��o na ra�z e ponta da se��o.

% beta [N.m] - vetor [1x2] que relaciona a rigidez � tors�o por unidade de
% largura (GJ/b) da se��o na ra�z e ponta da se��o com a rigidez a flex�o e
% a raz�o G/E da seguinte forma: K = beta * D * (G/E)

% b [m] - vetor [1x2] que indica a largura da se��o na ra�z e na ponta.

% ================================= BA ================================== %
% fBA [] - vetor [1x(k+1)] (onde k � o n�mero de se��es da longarina) que 
% indica a fra��o de corda que a longarina utiliza do ponto do BA at� o 
% ponto mais traseiro da longarina

% nx [] - escalar que indica o n�mero de pontos na dire��o x para a
% discretiza��o do perfil em cada se��o

% Cuidado: alguns par�metros globais n�o se aplicam a longarina, por�m
% mantenha-os para evitar possiveis erros. Os par�metros n�o utilizados 
% s�o: fc, S_fc, fh, l_long e d_long
% Centro�de � definido pela geometria da asa e n�o pelo criaLongarina


function [ L_w ] = ETT_criaLongarinaRelatorio( materiais )
%% INICIALIZA��O LONGARINA ASA %%
% ============================ INICIALIZA��O ============================ %
 L_w = struct;
 L_w(3).nome = 'Longarina  BA Asa Inferior';
 L_w(3).min_CS = 1.5;
 L_w(3).num_pontos = 200;
 L_w(3).is_long = true;
 L_w(3).is_tirante = false;
% % ============================ Asa Inferior ============================= %
 L_w(3).idx_malha = 1;
 L_w(3).param = "spline";
 L_w(3).config = {'BA','BA','BA','BA','BA','BA','BA','BA','BA','BA'};
 L_w(3).material = {materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X};
% 
%     % ======================== Par�metros de MEF ======================== %
     L_w(3).cc = { 3 , 1 , 2:2:6 }';
     L_w(3).idx_F = false;
     %conectar ate o pos 9
     L_w(3).connect = { 3 , 1 , @(L_w) L_w.posLong(2) , 'nearest';3 , 1 , @(L_w) L_w.posLong(3) , 'nearest';3 , 1 , @(L_w) L_w.posLong(4) , 'nearest';3 , 1 , @(L_w) L_w.posLong(5) , 'nearest';3 , 1 , @(L_w) L_w.posLong(6) , 'nearest';3 , 1 , @(L_w) L_w.posLong(7) , 'nearest';3 , 1 , @(L_w) L_w.posLong(8) , 'nearest';3 , 1 , @(L_w) L_w.posLong(9) , 'nearest';3 , 1 , @(L_w) L_w.posLong(10) , 'nearest';3 , 1 , @(L_w) L_w.posLong(end) , 'nearest'};
% 
%     % ======================== Par�metros Globais ======================= %
%         % ==================== Fra��o de Envergadura ==================== %
         L_w(3).fb = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
         L_w(3).fb0 = 0;
%         % ======================= Fra��o de Corda ======================= %
         L_w(3).fc = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
         L_w(3).S_fc = [ 0.2 , 0.2 ,0.2, 0.2 , 0.2, 0.2 ,0.2,0.2,0.2,0.2];
%         % ======================= Fra��o de Altura ====================== %
         L_w(3).fh = [0 0 0 0 0 0 0 0 0 0 0];
%         % ======================== Enflechamento ======================== %
         L_w(3).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
         L_w(3).d_long = [nan nan nan];
        
    % ======================== Par�metros BA ============================ %
        % ======================= Fra��o de Corda ======================= %
         L_w(3).fBA = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
%         % ================ N�mero de Discretiza��o em X ================= %
         L_w(3).nx = 100;
%         % ====================== N�mero de Voltas ======================= %
         L_w(3).n_voltas = [6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6];
%         % ====================== Fra��o com Voltas ====================== %
         L_w(3).f_voltas = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];

%% INICIALIZA��O LONGARINA SANDUICHE INFERIOR %%
% ======================= Asa Inferior - Sanduiche ======================= %
L_w(1).nome = 'Longarina Inferior Sanduiche H';
L_w(1).min_CS = 1.5;
L_w(1).num_pontos = 200;
L_w(1).is_long = true;
L_w(1).is_tirante = false;

L_w(1).idx_malha = 1;
L_w(1).config = {'sanduiche_I','sanduiche_I','sanduiche_I'};
L_w(1).material = {materiais.carbono,materiais.carbono,materiais.carbono};
L_w(1).material_c = {materiais.divinycellH45,materiais.divinycellH45,materiais.divinycellH45};

    % ======================== Par�metros de MEF ======================== %
    L_w(1).cc = { 1 , 1 , 1:6 }';
    L_w(1).idx_F = true;
%     L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
    L_w(1).connect = {};

    % ======================== Par�metros Globais ======================= %
        % ==================== Fra��o de Envergadura ==================== %
        L_w(1).fb = [0.217817500000000,0.450559900000000,1];
        L_w(1).fb0 = 0;
        % ======================= Fra��o de Corda ======================= %
        L_w(1).fc =[0.185325000000000,0.185325000000000,0.230655000000000,0.234960000000000];
        L_w(1).S_fc = [0.5 0.5 0.5];
        % ======================= Fra��o de Altura ====================== %
        L_w(1).fh = [0 0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(1).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(1).d_long = [nan nan nan];
        
    % ===================== Par�metros Sanduiche I ====================== %
        % ====================== Altura do N�cleo ======================= %
        L_w(1).H = [0.933418000000000,0.706060000000000,0.604189000000000,0.102538000000000];
        % ====================== N�mero de Camadas ====================== %
        L_w(1).n_lam = [1 1;1 1; 1 1];
        % ====================== Fra��o com Camadas ===================== %
        L_w(1).f_lam = [0.5 0.5 0.5];
        % =========================== Largura =========================== %
        L_w(1).B = [0.031795460000000,0.020698880000000,0.008531640000000,0.006241740000000];
%% INICIALIZA��O LONGARINAS - ASA SUPERIOR %%
% % ============================ INICIALIZA��O ============================ %
 L_w(4).nome = 'Longarina BA Asa Superior';
 L_w(4).min_CS = 1.5;
 L_w(4).num_pontos = 200;
 L_w(4).is_long = true;
 L_w(4).is_tirante = false;
% ============================ Asa Inferior ============================= %
 L_w(4).idx_malha = 2;
 L_w(4).param = "spline";
 L_w(4).config = {'BA','BA','BA','BA','BA','BA','BA','BA','BA','BA'};
 L_w(4).material = {materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X,materiais.balsa_X};
 
%     % ======================== Par�metros de MEF ======================== %
     L_w(4).cc = {4 , 1 , 2:2:6 }';
     L_w(4).idx_F = false;
     L_w(4).connect = { 4 , 2 , @(L_w) L_w.posLong(2), 'nearest';4 , 2 , @(L_w) L_w.posLong(3), 'nearest';4 , 2 , @(L_w) L_w.posLong(4), 'nearest';4 , 2 , @(L_w) L_w.posLong(5), 'nearest';4 , 2 , @(L_w) L_w.posLong(6), 'nearest';4 , 2 , @(L_w) L_w.posLong(7), 'nearest';4 , 2 , @(L_w) L_w.posLong(8), 'nearest';4 , 2 , @(L_w) L_w.posLong(9), 'nearest';4 , 2 , @(L_w) L_w.posLong(10), 'nearest';4 , 2 , @(L_w) L_w.posLong(end), 'nearest'};
% 
%     % ======================== Par�metros Globais ======================= %
%         % ==================== Fra��o de Envergadura ==================== %
         L_w(4).fb = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
         L_w(4).fb0 = 0;
%         % ======================= Fra��o de Corda ======================= %
         L_w(4).fc = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
         L_w(4).S_fc = [ 0.5 , 0.5, 0.5,0.5,0.5, 0.5, 0.5, 0.5 ,0.5, 0.5];
%         % ======================= Fra��o de Altura ====================== %
         L_w(4).fh = [0 0 0 0 0 0 0 0 0 0 0];
%         % ======================== Enflechamento ======================== %
         L_w(4).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
         L_w(4).d_long = [nan nan nan];
        
    % ======================== Par�metros BA ============================ %
%         % ======================= Fra��o de Corda ======================= %
         L_w(4).fBA = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; %ANALISAR O ERRO NO FBA
%         % ================ N�mero de Discretiza��o em X ================= %
         L_w(4).nx = 100;
%         % ====================== N�mero de Voltas ======================= %
         L_w(4).n_voltas = [6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6; 6 6];
%         % ====================== Fra��o com Voltas ====================== %
         L_w(4).f_voltas = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
%===============INICIALIZANDO SANDUICHE SUPERIOR===================%
L_w(2).nome = 'Longarina Superior Sanduiche H';
L_w(2).min_CS = 1.5;
L_w(2).num_pontos = 200;
L_w(2).is_long = true;
L_w(2).is_tirante = false;

L_w(2).idx_malha = 2;
L_w(2).config = {'sanduiche_I','sanduiche_I','sanduiche_I'};
L_w(2).material = {materiais.carbono,materiais.carbono,materiais.carbono};
L_w(2).material_c = {materiais.divinycellH45,materiais.divinycellH45,materiais.divinycellH45};

    % ======================== Par�metros de MEF ======================== %
    L_w(2).cc = { 2 , 1 , 1:6 }';
    L_w(2).idx_F = true;
%     L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
    L_w(2).connect = {};

    % ======================== Par�metros Globais ======================= %
        % ==================== Fra��o de Envergadura ==================== %
        L_w(2).fb = [0.112107499999966,0.509617599999973,1];
        L_w(2).fb0 = 0;
        % ======================= Fra��o de Corda ======================= %
        L_w(2).fc =[0.214362500000000,0.214362500000000,0.200470000000000,0.184345000000000];
        L_w(2).S_fc = [0.5 0.5 0.5];
        % ======================= Fra��o de Altura ====================== %
        L_w(2).fh = [0 0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(2).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(2).d_long = [nan nan nan];
        
    % ===================== Par�metros Sanduiche I ====================== %
        % ====================== Altura do N�cleo ======================= %
        L_w(2).H = [0.80569,0.805690000000000,0.755272000000000,0.300088000000000];
        % ====================== N�mero de Camadas ====================== %
        L_w(2).n_lam = [1 1;1 1; 1 1];
        % ====================== Fra��o com Camadas ===================== %
        L_w(2).f_lam = [0.5 0.5 0.5];
        % =========================== Largura =========================== %
        L_w(2).B = [0.02470120000000,0.018688120000000,0.007726860000000,0.006000000000000];
%% INICIALIZA��O TIRANTE POSTERIOR %%
% ============================ INICIALIZA��O ============================ %
% L_w(3).nome = 'Tirante Posterior';
% L_w(3).min_CS = 1.5;
% L_w(3).num_pontos = 20;
% L_w(3).is_long = false;
% L_w(3).is_tirante = true;
% % ============================ Asa Inferior ============================= %
% L_w(3).idx_malha = [1 2];
% L_w(3).config = {'sanduiche'};
% L_w(3).material = {materiais.carbono};
% 
%     % ======================== Par�metros de MEF ======================== %
%     L_w(3).cc = {};
%     L_w(3).idx_F = false;
%     L_w(3).connect = { 3 , 1 , @(L_w) L_w.posLong(1) , 'nearest' ; 3 , 2 , @(L_w) L_w.posLong(end) , 'nearest' };
% 
%     % ======================== Par�metros Tirante ======================= %
%         % ======================= Fra��o de Corda ======================= %
%         L_w(3).fc = [0.20 0.35];
%         % ==================== Fra��o de Envergadura ==================== %
%         L_w(3).fb = [0 0];
%         
%     % ====================== Par�metros Sandu�che ======================= %
%         % ======================= Material N�cleo ======================= %
%         L_w(3).material_c = {materiais.divinycellH60};
%         % ====================== In�rcia a Flex�o ======================= %
%         L_w(3).D_I = {[7000 7000 7000]};
%         % ============= Fator Beta (define in�rcia a tors�o) ============ %
%         L_w(3).beta = {[1 1 0.5]};
%         % ====================== Largura da Se��o ======================= %
%         L_w(3).b = {[0.08 0.08 0.05]};


%% INICIALIZA��O TIRANTE ANTERIOR %%
% ============================ INICIALIZA��O ============================ %
% L_w(7).nome = 'Fios';
% L_w(7).min_CS = 1.5;
% L_w(7).num_pontos = 10;
% L_w(7).is_long = false;
% L_w(7).is_tirante = true;
% % ============================ Asa Inferior ============================= %
% L_w(7).idx_malha = [1 2];
% L_w(7).config = {'fio'};
% L_w(7).material = {materiais.kevlar};
% 
%     %======================== Par�metros de MEF ======================== %
%     L_w(7).cc = {};
%     L_w(7).idx_F = false;
%     L_w(7).connect = { 7 , 1 , @(L_w) L_w.posLong(1) , 'nearest' ; 7 , 2 , @(L_w) L_w.posLong(end) , 'nearest' };
%  
%     %======================== Par�metros Tirante ======================= %
%         %======================= Fra��o de Corda ======================= %
%         L_w(7).fc = [0.3 0.3];
%        % ==================== Fra��o de Envergadura ==================== %
%         L_w(7).fb = [0.8 0.8];
%         
%    % ====================== Par�metros Cicular ======================= %
%     L_w(7).D = [0.001 0.001];    
% % ======================= Empenagem Horizontal - Sanduiche ======================= %
%L_w(5).nome = 'Longarina Empenagem Horizontal';
%L_w(5).min_CS = 1.5;
%L_w(5).num_pontos = 100;
%L_w(5).is_long = true;
%L_w(5).is_tirante = false;

%L_w(5).idx_malha = 3;
%L_w(5).config = {'sanduiche_H','sanduiche_H','sanduiche_H'};
%L_w(5).material = {materiais.carbono,materiais.carbono,materiais.carbono};
%L_w(5).material_c = {materiais.divinycellH60,materiais.divinycellH60,materiais.divinycellH60};

    % ======================== Par�metros de MEF ======================== %
    %L_w(5).cc = { 5 , 1 , 1:6 }';
    %L_w(5).idx_F = true;
%     L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
    %L_w(5).connect = {};

    % ======================== Par�metros Globais ======================= %
        % ==================== Fra��o de Envergadura ==================== %
        %L_w(5).fb = [0.2 0.6 1];
        %L_w(5).fb0 = 0;
        % ======================= Fra��o de Corda ======================= %
        %L_w(5).fc =[0.7 0.7 0.7 0.7];
        %L_w(5).S_fc = [0.5 0.5 0.5];
        % ======================= Fra��o de Altura ====================== %
        %L_w(5).fh = [0 0 0 0];
        % ======================== Enflechamento ======================== %
        %L_w(5).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        %L_w(5).d_long = [nan nan nan];
        
    % ===================== Par�metros Sanduiche I ====================== %
        % ====================== Altura do N�cleo ======================= %
        %L_w(5).H = [1 1 1 1];
        % ====================== N�mero de Camadas ====================== %
        %L_w(5).n_lam = [1 1;1 1; 1 1];
        % ====================== Fra��o com Camadas ===================== %
        %L_w(5).f_lam = [0.5 0.5 0.5];
        % =========================== Largura =========================== %
        %L_w(5).B = [0.01 0.01 0.01 0.01];
% %         
% % =============== Longarina BA Empenagem Horizontal ================= %
 %L_w(6).nome = 'Longarina  BA Empenagem Horizontal';
 %L_w(6).min_CS = 1.5;
 %L_w(6).num_pontos = 100;
 %L_w(6).is_long = true;
 %L_w(6).is_tirante = false;
% % ============================ Asa Inferior ============================= %
 %L_w(6).idx_malha = 3;
 %L_w(6).param = "discreto";
 %L_w(6).config = {'BA','BA'};
 %L_w(6).material = {materiais.balsa_X,materiais.balsa_X};

    % ======================== Par�metros de MEF ======================== %
     %L_w(6).cc = { 6 , 1 , 2:2:6 }';
     %L_w(6).idx_F = true;
     %L_w(6).connect = { 6 , 5 , @(L_w) L_w.posLong(2) , 'nearest';6 , 5 , @(L_w) L_w.posLong(end) , 'nearest'};

    % ======================== Par�metros Globais ======================= %
        % ==================== Fra��o de Envergadura ==================== %
         %L_w(6).fb = [0.4912 1]; %COLOCAR O FB DA MALHA DA ASA (N�o � parametro da longarina)
         %L_w(6).fb0 = 0;
        % ======================= Fra��o de Corda ======================= %
         %L_w(6).fc = [0.1 0.1 0.1];
         %L_w(6).S_fc = [ 0.2 , 0.2];
        % ======================= Fra��o de Altura ====================== %
         %L_w(6).fh = [0 0 0];
        % =======================3 , 1 , @(L_w) L_w.posLong(5) , 'nearest'= Enflechamento ======================== %
         %L_w(6).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
         %L_w(6).d_long = [nan nan nan];
        
    % ======================== Par�metros BA ============================ %
       % ======================= Fra��o de Corda ======================= %
         %L_w(6).fBA = [0.2 0.2 0.2];
        % ================ N�mero de Discretiza��o em X ================= %
         %L_w(6).nx = 1000;
        % ====================== N�mero de Voltas ======================= %
         %L_w(6).n_voltas = [6 6;6 6];
        % ====================== Fra��o com Voltas ====================== %
         %L_w(6).f_voltas = [0.5 0.5];
% %% INICIALIZA��O TUBINHO POSTERIOR %%
% % ============================ INICIALIZA��O ============================ %
%  L_w(5).nome = 'Tubinho Asa Inferior';
%  L_w(5).min_CS = 1.5;
%  L_w(5).num_pontos = 100;
%  L_w(5).is_long = true;
%  L_w(5).is_tirante = false;
% % ============================ Asa Inferior ============================= %
%  L_w(5).idx_malha = 3;
%  L_w(5).config = {'tubinho'};
%  L_w(5).material = {materiais.carbono};
% 
%     % ======================== Par�metros de MEF ======================== %
%      L_w(5).cc = {5 , 1, 1:6}';
%      L_w(5).idx_F = true;
%      L_w(5).connect = {};
% 
%     % ======================== Par�metros Globais ======================= %
%         % ==================== Fra��o de Envergadura ==================== %
%          L_w(5).fb = [1];
%          L_w(5).fb0 = 0;
%         % ======================= Fra��o de Corda ======================= %
%          L_w(5).fc = [ 0.11 0.11 ];
%          L_w(5).S_fc = [ nan ];
%         % ======================= Fra��o de Altura ====================== %
%          L_w(5).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%          L_w(5).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%          L_w(5).d_long = [nan nan nan];
%          % ======================= Par�metros Tubinho ======================== %
%        % ====================== Di�metro Interno ======================= %
%          L_w(5).D_interno = [2] * 1e-3;
%         % ====================== Di�metro Externo ======================= %
%          L_w(5).D_externo = [3] * 1e-3;
%          
%  %% ========== SEGUNDO TUBINHO =============== %%       
%  L_w(6).nome = 'Tubinho Asa Inferior';
%  L_w(6).min_CS = 1.5;
%  L_w(6).num_pontos = 100;
%  L_w(6).is_long = true;
%  L_w(6).is_tirante = false;
% % ============================ Asa Inferior ============================= %
%  L_w(6).idx_malha = 3;
%  L_w(6).config = {'tubinho'};
%  L_w(6).material = {materiais.carbono};
% 
%     % ======================== Par�metros de MEF ======================== %
%      L_w(6).cc = {6 , 1, 1:6}';
%      L_w(6).idx_F = true;
%      L_w(6).connect = { 6 , 5 , @(L_w) L_w.posLong(end) , 'nearest'};
% 
%     % ======================== Par�metros Globais ======================= %
%         % ==================== Fra��o de Envergadura ==================== %
%          L_w(6).fb = [1];
%          L_w(6).fb0 = 0;
%         % ======================= Fra��o de Corda ======================= %
%          L_w(6).fc = [ 0.31 0.31 ];
%          L_w(6).S_fc = [ nan ];
%         % ======================= Fra��o de Altura ====================== %
%          L_w(6).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%          L_w(6).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%          L_w(6).d_long = [nan nan nan];
%    % ======================= Par�metros Tubinho ======================== %
%        % ====================== Di�metro Interno ======================= %
%          L_w(6).D_interno = [2] * 1e-3;
%         % ====================== Di�metro Externo ======================= %
%          L_w(6).D_externo = [3] * 1e-3;

% %% INICIALIZA��O TUBINHO POSTERIOR %%
% % ============================ INICIALIZA��O ============================ %
% L_w(7).nome = 'Tubinho Asa Superior';
% L_w(7).min_CS = 1.5;
% L_w(7).num_pontos = 10;
% L_w(7).is_long = true;
% L_w(7).is_tirante = false;
% % ============================ Asa Inferior ============================= %
% L_w(7).idx_malha = 2;
% L_w(7).config = {'tubinho'};
% L_w(7).material = {materiais.carbono};
% 
%     % ======================== Par�metros de MEF ======================== %
%     L_w(7).cc = {}';
%     L_w(7).idx_F = true;
%     L_w(7).connect = { 7 , 3 , @(L_w) L_w.posLong(1) , 'nearest' };
% 
%     % ======================== Par�metros Globais ======================= %
%         % ==================== Fra��o de Envergadura ==================== %
%         L_w(7).fb = [1];
%         L_w(7).fb0 = 1 - prod(1 - L_w(3).fb);
%         % ======================= Fra��o de Corda ======================= %
%         L_w(7).fc = [ 0.3 0.32 ];
%         L_w(7).S_fc = [ nan ];
%         % ======================= Fra��o de Altura ====================== %
%         L_w(7).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%         L_w(7).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%         L_w(7).d_long = [nan nan nan];
%         
%     % ======================= Par�metros Tubinho ======================== %
%         % ====================== Di�metro Interno ======================= %
%         L_w(7).D_interno = [6] * 1e-3;
%         % ====================== Di�metro Externo ======================= %
%         L_w(7).D_externo = [8] * 1e-3;

end

