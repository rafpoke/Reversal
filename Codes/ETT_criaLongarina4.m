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


function [ L_w ] = ETT_criaLongarina4( materiais )
%% INICIALIZA��O LONGARINA ASA %%
% ============================ INICIALIZA��O ============================ %
L_w = struct;
L_w(1).nome = 'Longarina Asa Inferior';
L_w(1).min_CS = 1.5;
L_w(1).num_pontos = 200;
L_w(1).is_long = true;
L_w(1).is_tirante = false;
L_w(2).is_trelica = false;
% ============================ Asa Inferior ============================= %
L_w(1).idx_malha = 1;
L_w(1).config = {'BA','BA'};
L_w(1).material = {materiais.carbono,materiais.carbono};

    % ======================== Par�metros de MEF ======================== %
    L_w(1).cc = { 1 , 1 , 1:6 }';
    L_w(1).idx_F = true;
    L_w(1).connect = {};

    % ======================== Par�metros Globais ======================= %
        % ==================== Fra��o de Envergadura ==================== %
        L_w(1).fb = [0.5 1];
        L_w(1).fb0 = 0;
        % ======================= Fra��o de Corda ======================= %
        L_w(1).fc = [0.4 0.4 0.2];
        L_w(1).S_fc = [ 0.5 , 0.2 ];
        % ======================= Fra��o de Altura ====================== %
        L_w(1).fh = [0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(1).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(1).d_long = [nan nan nan];
        
    % ======================== Par�metros BA ============================ %
        % ======================= Fra��o de Corda ======================= %
        L_w(1).fBA = [0.2 0.1 0.05];
        % ================ N�mero de Discretiza��o em X ================= %
        L_w(1).nx = 100;
        % ====================== N�mero de Voltas ======================= %
        L_w(1).n_voltas = [1 1;1 1];
        % ====================== Fra��o com Voltas ====================== %
        L_w(1).f_voltas = [0.5 0.5];

%% INICIALIZA��O LONGARINA AUXILIAR %%
% ======================= Asa Inferior - Auxiliar ======================= %
L_w(2).nome = 'Longarina Auxiliar';
L_w(2).min_CS = 1.5;
L_w(2).num_pontos = 100;
L_w(2).is_long = true;
L_w(2).is_tirante = false;
L_w(2).is_trelica = false;

L_w(2).idx_malha = 1;
L_w(2).config = {'circular'};
L_w(2).material = {materiais.carbono};

    % ======================== Par�metros de MEF ======================== %
    L_w(2).cc = { 2 , 1 , 1:6 }';
    L_w(2).idx_F = true;
    L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
%     L_w(2).connect = {};

    % ======================== Par�metros Globais ======================= %
        % ==================== Fra��o de Envergadura ==================== %
        L_w(2).fb = [0.22];
        L_w(2).fb0 = 0;
        % ======================= Fra��o de Corda ======================= %
        L_w(2).fc = [0.2 0.2];
        L_w(2).S_fc = 0.5;
        % ======================= Fra��o de Altura ====================== %
        L_w(2).fh = [0 0];
        % ======================== Enflechamento ======================== %
        L_w(2).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(2).d_long = [nan nan nan];
        
    % ======================= Par�metros Circular ======================= %
        % ========================== Di�metro =========================== %
        L_w(2).D = [9/8] * 25.4e-3;
        % ====================== N�mero de Voltas ======================= %
        L_w(2).n_voltas = [3 2];
        % ====================== Fra��o com Voltas ====================== %
        L_w(2).f_voltas = [0.5];
        
%% INICIALIZA��O LONGARINAS - ASA SUPERIOR %%
% ============================ INICIALIZA��O ============================ %
L_w(3).nome = 'Longarina Asa Superior';
L_w(3).min_CS = 1.5;
L_w(3).num_pontos = 200;
L_w(3).is_long = true;
L_w(3).is_tirante = false;
L_w(3).is_trelica = false;
% ============================ Asa Inferior ============================= %
L_w(3).idx_malha = 2;
L_w(3).config = {'BA','BA'};
L_w(3).material = {materiais.carbono,materiais.carbono};

    % ======================== Par�metros de MEF ======================== %
    L_w(3).cc = { 1 , 1 , 1:6 }';
    L_w(3).idx_F = true;
    L_w(3).connect = {};

    % ======================== Par�metros Globais ======================= %
        % ==================== Fra��o de Envergadura ==================== %
        L_w(3).fb = [0.5 1];
        L_w(3).fb0 = 0;
        % ======================= Fra��o de Corda ======================= %
        L_w(3).fc = [0.4 0.4 0.2];
        L_w(3).S_fc = [ 0.5 , 0.2 ];
        % ======================= Fra��o de Altura ====================== %
        L_w(3).fh = [0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(3).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(3).d_long = [nan nan nan];
        
    % ======================== Par�metros BA ============================ %
        % ======================= Fra��o de Corda ======================= %
        L_w(3).fBA = [0.3 0.1 0.05];
        % ================ N�mero de Discretiza��o em X ================= %
        L_w(3).nx = 100;
        % ====================== N�mero de Voltas ======================= %
        L_w(3).n_voltas = [1 1;1 1];
        % ====================== Fra��o com Voltas ====================== %
        L_w(3).f_voltas = [0.5 0.5];

        
%% INICIALIZA��O TIRANTE POSTERIOR %%
% ============================ INICIALIZA��O ============================ %
% L_w(4).nome = 'Tirante Posterior';
% L_w(4).min_CS = 1.5;
% L_w(4).num_pontos = 25;
% L_w(4).is_long = false;
% L_w(4).is_tirante = true;
% L_w(4).is_trelica = true;
% ============================ Asa Inferior ============================= %
% L_w(4).idx_malha = [1 2];
% L_w(4).config = {'quadrado'};
% L_w(4).material = {materiais.carbono,materiais.carbono};

    % ======================== Par�metros de MEF ======================== %
   % L_w(4).cc = {};
   % L_w(4).idx_F = false;
   % L_w(4).connect = { 4 , 2 , @(L_w) L_w.posLong(1) , 'nearest' ; 4 , 3 , @(L_w) L_w.posLong(end) , 'nearest' };

    % ======================== Par�metros Tirante ======================= %
        % ======================= Fra��o de Corda ======================= %
       % L_w(4).fc = [0.15 0.15];
        % ==================== Fra��o de Envergadura ==================== %
       % L_w(4).fb = [0 0];
        % ======================= N�mero de se��es ====================== %
       % L_w(4).num_T = 5;
    % ======================== Par�metros Treli�a ======================= %
        % =========================== Largura =========================== %
       % L_w(4).b = [0.08 0.08;0.08 0.08];
        % ==================== Curvatura da Largura ===================== %
       % L_w(4).Sb = [0.8642 0.8642];
        % ======================= Comprimento em X ====================== %
       % L_w(4).h = [0.08 0.04 0.08];
        % =================== Curvatura do Comprimento ================== %
       % L_w(4).Sh = 0.92;


%% INICIALIZA��O TIRANTE ANTERIOR %%
% ============================ INICIALIZA��O ============================ %
% L_w(5).nome = 'Tirante Posterior';
% L_w(5).min_CS = 1.5;
% L_w(5).num_pontos = 6;
% L_w(5).is_long = false;
% L_w(5).is_tirante = true;
% % ============================ Asa Inferior ============================= %
% L_w(5).idx_malha = [1 2];
% L_w(5).config = {'sanduiche'};
% L_w(5).material = {materiais.carbono};
% 
%     % ======================== Par�metros de MEF ======================== %
%     L_w(5).cc = {};
%     L_w(5).idx_F = false;
%     L_w(5).connect = { 5 , 1 , @(L_w) L_w.posLong(1) , 'nearest' ; 5 , 3 , @(L_w) L_w.posLong(end) , 'nearest' };
% 
%     % ======================== Par�metros Tirante ======================= %
%         % ======================= Fra��o de Corda ======================= %
%         L_w(5).fc = [0.45 0.3];
%         % ==================== Fra��o de Envergadura ==================== %
%         L_w(5).fb = [0 0];
%         
%     % ====================== Par�metros Sandu�che ======================= %
%         % ======================= Material N�cleo ======================= %
%         L_w(5).material_c = {materiais.divinycell};
%         % ====================== In�rcia a Flex�o ======================= %
%         L_w(5).D_I = {[600 600 0.5]};
%         % ============= Fator Beta (define in�rcia a tors�o) ============ %
%         L_w(5).beta = {[1 1 0.5];[1 1 0.5]};
%         % ====================== Largura da Se��o ======================= %
%         L_w(5).b = {[0.05 0.05 0.5]};
%         
% %% INICIALIZA��O TUBINHO POSTERIOR %%
% % ============================ INICIALIZA��O ============================ %
% L_w(6).nome = 'Tubinho Asa Inferior';
% L_w(6).min_CS = 1.5;
% L_w(6).num_pontos = 10;
% L_w(6).is_long = true;
% L_w(6).is_tirante = false;
% % ============================ Asa Inferior ============================= %
% L_w(6).idx_malha = 1;
% L_w(6).config = {'tubinho'};
% L_w(6).material = {materiais.carbono};
% 
%     % ======================== Par�metros de MEF ======================== %
%     L_w(6).cc = {}';
%     L_w(6).idx_F = true;
%     L_w(6).connect = { 6 , 1 , @(L_w) L_w.posLong(1) , 'nearest' };
% 
%     % ======================== Par�metros Globais ======================= %
%         % ==================== Fra��o de Envergadura ==================== %
%         L_w(6).fb = [1];
%         L_w(6).fb0 = 1 - prod(1 - L_w(1).fb);
%         % ======================= Fra��o de Corda ======================= %
%         L_w(6).fc = [ 0.355 0.32 ];
%         L_w(6).S_fc = [ nan ];
%         % ======================= Fra��o de Altura ====================== %
%         L_w(6).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%         L_w(6).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%         L_w(6).d_long = [nan nan nan];
%         
%     % ======================= Par�metros Tubinho ======================== %
%         % ====================== Di�metro Interno ======================= %
%         L_w(6).D_interno = [2] * 1e-3;
%         % ====================== Di�metro Externo ======================= %
%         L_w(6).D_externo = [4] * 1e-3;
% 
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

