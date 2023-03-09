%% Função ETT_criaLongarina %%
% ============================== DESCRIÇÃO ============================== %
% Essa função parametriza completamente a geometria de longarinas com
% múltiplas seções. Essa parametrização é útil no caso de integração com
% uma rotina de otimização. Os parâmetros podem ser de dois tipos: globais
% ou específicos. Parâmetros globais são utilizados para descrever qualquer
% tipo de longarina (independentemente da geometria da seção). Parâmetros
% específicos descrevem uma geometria específica (como circular, retangular
% etc). Para auxiliar o entendimento da função é recomendado alterar
% parâmetros e utilizar a função plotLongarina() para melhor visualização.

% =============================== OUTPUTS =============================== %
% - Longarina - struct [1xn] onde n é o número de longarinas parametrizadas
% da aeronave.

% =========================== PARAMETRIZAÇÃO ============================ %

% =============================== GLOBAL ================================ %
% - config - cell array [1xk] (onde k é o número de seções da longarina) de
% strings. As string simbolizam o formato da seção da longarina. Formatos
% possíveis: 'circular' (programado), 'retangular' (não-programado), 'I'
% (não-programado), 'C' (não-programado), 'BA' (não-programado), 'tubo'
% (não-programado).
% Ex: Longarina(1).config = {'circular','retangular'}

% - material - cell array [1xk] (onde k é o número de seções da longarina)
% de structs. As structs indicam o material do qual é feita a seção de
% longarina. Materiais possíveis: carbono (programado), fibra_vidro
% (não-programado), dyneema (não-programado), kevlar_carbono 
% (não-programado), balsa (não-programado), plywood (não-programado).
% Ex: Longarina(3).material = {materiais.carbono,materiais.balsa}

% - cc - cell array [nx2] (onde n é o número de condições de contorno da
% longarina). Cada linha indica o número do nó em que se aplica a condição
% de contorno e os graus de liberdade restringidos de 1 a 6 (u v w tx ty
% tz).

% - fb - vetor [1x(k-1)] (onde k é o número de seções da longarina) da
% fração da envergadura restante ocupada pela seção da longarina. Por
% exemplo em uma asa de semi_envergadura b=1 e fb = [0.5 0.5] teríamos uma
% longarina com três seções de envergaduras [0.5 0.25 0.25]; no caso de fb
% = [0.4 0.3] teríamos envergaduras [0.4 0.18 0.42].
% Ex: Longarina(1).fb = [0.4] ou Longarina(1).fb = [] (no caso de uma seção
% única de longarina)

% - fc - vetor [1x(k+1)] (onde k é o número de seções da longarina) que 
% indicam a fração de corda em que estão posicionados os pontos de
% troca de seção da longarina.
% Ex: Longarina(2).fc = [0.5 0.4 0.25] (nesse caso a seção 1 teria
% extremidades em 50% de corda e 40% de corda, enquanto a seção 2 teria
% extremidades em 40% de corda e 25% de corda)

% - fh - vetor [1x(k+1)] (onde  k é o número de seções da longarina) que
% indica o quão afastado em z do centro do perfil está o centro da seção
% da longarina. Esse valor é uma fração da distância em z da parte superior
% da asa à parte inferior que pode variar entre [-0.5,+0.5].
% Ex: Longarina(1).fh = [0 -0.05 0.1]

% ============================== CIRCULAR =============================== %
% - fr - vetor [1xk] (onde  k é o número de seções da longarina) que indica
% a fração do raio máximo da seção que será utlizado. Por exemplo, se a asa
% comportaria uma seção j de longarina circular com raio 15mm e fr(j) = 0.5
% então a seção j da longarina terá raio de 7.5mm.
% Ex: Longarina(2).fr = [0.7 0.6]

% - f_voltas - vetor [1xk] (onde  k é o número de seções da longarina) que
% indica a fração de envergadura da seção em que ocorrerá uma mudança de
% número de voltas de laminado. Por exemplo, se a envergadura da seção j é
% b=0.6 e f_voltas(j)=0.4 então ocorrerá uma troca de número de voltas em
% y=0.24.
% Ex: Longarina(2).f_voltas = [0.2 nan 0.8] (no caso das seções [1,3] serem
% circulares e a seção [2] não)

% - n_voltas - matriz [kx2] (onde  k é o número de seções da longarina) em
% que cada linha indica o número de voltas de laminado de cada seção de
% longarina.
% Ex: Longarina(3).n_voltas = [3 2;nan nan;2 2] (nesse caso na primeira
% seção a longarina inicialmente tem 3 voltas e depois tem 2. A segunda
% seção não é circular e a terceira tem 2 voltas ao longo de toda a
% envergadura)

% ============================= RETANGULAR ============================== %
% phi [deg] - vetor [1x(k+1)] (onde k é o número de seções da longarina)
% que indica a relação entre a alma e a base (phi=atand(alma/base)) em cada
% troca de seção da longarina.
% Ex: Longarina(1) = [45 45 45] (nesse caso a longarina tem 2 seções, todas
% com alma igual à base).

% e_alma [m] - vetor [1xk] (onde k é o número de seções da longarina) que
% indica a espessura da alma da longarina retangular em cada seção. A alma 
% é a lateral da longarina.

% e_base [m] - vetor [1xk] (onde k é o número de seções da longarina) que
% indica a espessura da base da longarina retangular em cada seção. A base 
% é a parte superior e inferior da longarina.

% f_dim - vetor [1x(k+1)] (onde k é o número de seções da longarina) que 
% indica a fração da dimensão máxima disponível em uma troca de seção que
% a longarina ocupará.
% Ex: Longarina(1).f_dim = [1 0.9] (nesse caso a longarina possui apenas
% uma seção. Se na raíz a dimensão máxima disponível para a alma fosse 40mm
% e na ponta fosse 10mm, as dimensões reais da alma seriam 40mm na raíz e
% 9mm na ponta.

% ============================== SANDUICHE ============================== %
% D_I [N.m] - vetor [1x2] que indica a rigidez à flexão por unidade de
% largura (EI/b) da seção na raíz e ponta da seção.

% beta [N.m] - vetor [1x2] que relaciona a rigidez à torsão por unidade de
% largura (GJ/b) da seção na raíz e ponta da seção com a rigidez a flexão e
% a razão G/E da seguinte forma: K = beta * D * (G/E)

% b [m] - vetor [1x2] que indica a largura da seção na raíz e na ponta.

% ================================= BA ================================== %
% fBA [] - vetor [1x(k+1)] (onde k é o número de seções da longarina) que 
% indica a fração de corda que a longarina utiliza do ponto do BA até o 
% ponto mais traseiro da longarina

% nx [] - escalar que indica o número de pontos na direção x para a
% discretização do perfil em cada seção

% Cuidado: alguns parâmetros globais não se aplicam a longarina, porém
% mantenha-os para evitar possiveis erros. Os parâmetros não utilizados 
% são: fc, S_fc, fh, l_long e d_long
% Centroíde é definido pela geometria da asa e não pelo criaLongarina


function [ L_w ] = ETT_criaLongarinaRelatorio( materiais )
%% INICIALIZAÇÃO LONGARINA ASA %%
% ============================ INICIALIZAÇÃO ============================ %
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
%     % ======================== Parâmetros de MEF ======================== %
     L_w(3).cc = { 3 , 1 , 2:2:6 }';
     L_w(3).idx_F = false;
     %conectar ate o pos 9
     L_w(3).connect = { 3 , 1 , @(L_w) L_w.posLong(2) , 'nearest';3 , 1 , @(L_w) L_w.posLong(3) , 'nearest';3 , 1 , @(L_w) L_w.posLong(4) , 'nearest';3 , 1 , @(L_w) L_w.posLong(5) , 'nearest';3 , 1 , @(L_w) L_w.posLong(6) , 'nearest';3 , 1 , @(L_w) L_w.posLong(7) , 'nearest';3 , 1 , @(L_w) L_w.posLong(8) , 'nearest';3 , 1 , @(L_w) L_w.posLong(9) , 'nearest';3 , 1 , @(L_w) L_w.posLong(10) , 'nearest';3 , 1 , @(L_w) L_w.posLong(end) , 'nearest'};
% 
%     % ======================== Parâmetros Globais ======================= %
%         % ==================== Fração de Envergadura ==================== %
         L_w(3).fb = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
         L_w(3).fb0 = 0;
%         % ======================= Fração de Corda ======================= %
         L_w(3).fc = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
         L_w(3).S_fc = [ 0.2 , 0.2 ,0.2, 0.2 , 0.2, 0.2 ,0.2,0.2,0.2,0.2];
%         % ======================= Fração de Altura ====================== %
         L_w(3).fh = [0 0 0 0 0 0 0 0 0 0 0];
%         % ======================== Enflechamento ======================== %
         L_w(3).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
         L_w(3).d_long = [nan nan nan];
        
    % ======================== Parâmetros BA ============================ %
        % ======================= Fração de Corda ======================= %
         L_w(3).fBA = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
%         % ================ Número de Discretização em X ================= %
         L_w(3).nx = 100;
%         % ====================== Número de Voltas ======================= %
         L_w(3).n_voltas = [6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6];
%         % ====================== Fração com Voltas ====================== %
         L_w(3).f_voltas = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];

%% INICIALIZAÇÃO LONGARINA SANDUICHE INFERIOR %%
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

    % ======================== Parâmetros de MEF ======================== %
    L_w(1).cc = { 1 , 1 , 1:6 }';
    L_w(1).idx_F = true;
%     L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
    L_w(1).connect = {};

    % ======================== Parâmetros Globais ======================= %
        % ==================== Fração de Envergadura ==================== %
        L_w(1).fb = [0.217817500000000,0.450559900000000,1];
        L_w(1).fb0 = 0;
        % ======================= Fração de Corda ======================= %
        L_w(1).fc =[0.185325000000000,0.185325000000000,0.230655000000000,0.234960000000000];
        L_w(1).S_fc = [0.5 0.5 0.5];
        % ======================= Fração de Altura ====================== %
        L_w(1).fh = [0 0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(1).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(1).d_long = [nan nan nan];
        
    % ===================== Parâmetros Sanduiche I ====================== %
        % ====================== Altura do Núcleo ======================= %
        L_w(1).H = [0.933418000000000,0.706060000000000,0.604189000000000,0.102538000000000];
        % ====================== Número de Camadas ====================== %
        L_w(1).n_lam = [1 1;1 1; 1 1];
        % ====================== Fração com Camadas ===================== %
        L_w(1).f_lam = [0.5 0.5 0.5];
        % =========================== Largura =========================== %
        L_w(1).B = [0.031795460000000,0.020698880000000,0.008531640000000,0.006241740000000];
%% INICIALIZAÇÃO LONGARINAS - ASA SUPERIOR %%
% % ============================ INICIALIZAÇÃO ============================ %
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
 
%     % ======================== Parâmetros de MEF ======================== %
     L_w(4).cc = {4 , 1 , 2:2:6 }';
     L_w(4).idx_F = false;
     L_w(4).connect = { 4 , 2 , @(L_w) L_w.posLong(2), 'nearest';4 , 2 , @(L_w) L_w.posLong(3), 'nearest';4 , 2 , @(L_w) L_w.posLong(4), 'nearest';4 , 2 , @(L_w) L_w.posLong(5), 'nearest';4 , 2 , @(L_w) L_w.posLong(6), 'nearest';4 , 2 , @(L_w) L_w.posLong(7), 'nearest';4 , 2 , @(L_w) L_w.posLong(8), 'nearest';4 , 2 , @(L_w) L_w.posLong(9), 'nearest';4 , 2 , @(L_w) L_w.posLong(10), 'nearest';4 , 2 , @(L_w) L_w.posLong(end), 'nearest'};
% 
%     % ======================== Parâmetros Globais ======================= %
%         % ==================== Fração de Envergadura ==================== %
         L_w(4).fb = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
         L_w(4).fb0 = 0;
%         % ======================= Fração de Corda ======================= %
         L_w(4).fc = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
         L_w(4).S_fc = [ 0.5 , 0.5, 0.5,0.5,0.5, 0.5, 0.5, 0.5 ,0.5, 0.5];
%         % ======================= Fração de Altura ====================== %
         L_w(4).fh = [0 0 0 0 0 0 0 0 0 0 0];
%         % ======================== Enflechamento ======================== %
         L_w(4).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
         L_w(4).d_long = [nan nan nan];
        
    % ======================== Parâmetros BA ============================ %
%         % ======================= Fração de Corda ======================= %
         L_w(4).fBA = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; %ANALISAR O ERRO NO FBA
%         % ================ Número de Discretização em X ================= %
         L_w(4).nx = 100;
%         % ====================== Número de Voltas ======================= %
         L_w(4).n_voltas = [6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6;6 6; 6 6];
%         % ====================== Fração com Voltas ====================== %
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

    % ======================== Parâmetros de MEF ======================== %
    L_w(2).cc = { 2 , 1 , 1:6 }';
    L_w(2).idx_F = true;
%     L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
    L_w(2).connect = {};

    % ======================== Parâmetros Globais ======================= %
        % ==================== Fração de Envergadura ==================== %
        L_w(2).fb = [0.112107499999966,0.509617599999973,1];
        L_w(2).fb0 = 0;
        % ======================= Fração de Corda ======================= %
        L_w(2).fc =[0.214362500000000,0.214362500000000,0.200470000000000,0.184345000000000];
        L_w(2).S_fc = [0.5 0.5 0.5];
        % ======================= Fração de Altura ====================== %
        L_w(2).fh = [0 0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(2).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(2).d_long = [nan nan nan];
        
    % ===================== Parâmetros Sanduiche I ====================== %
        % ====================== Altura do Núcleo ======================= %
        L_w(2).H = [0.80569,0.805690000000000,0.755272000000000,0.300088000000000];
        % ====================== Número de Camadas ====================== %
        L_w(2).n_lam = [1 1;1 1; 1 1];
        % ====================== Fração com Camadas ===================== %
        L_w(2).f_lam = [0.5 0.5 0.5];
        % =========================== Largura =========================== %
        L_w(2).B = [0.02470120000000,0.018688120000000,0.007726860000000,0.006000000000000];
%% INICIALIZAÇÃO TIRANTE POSTERIOR %%
% ============================ INICIALIZAÇÃO ============================ %
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
%     % ======================== Parâmetros de MEF ======================== %
%     L_w(3).cc = {};
%     L_w(3).idx_F = false;
%     L_w(3).connect = { 3 , 1 , @(L_w) L_w.posLong(1) , 'nearest' ; 3 , 2 , @(L_w) L_w.posLong(end) , 'nearest' };
% 
%     % ======================== Parâmetros Tirante ======================= %
%         % ======================= Fração de Corda ======================= %
%         L_w(3).fc = [0.20 0.35];
%         % ==================== Fração de Envergadura ==================== %
%         L_w(3).fb = [0 0];
%         
%     % ====================== Parâmetros Sanduíche ======================= %
%         % ======================= Material Núcleo ======================= %
%         L_w(3).material_c = {materiais.divinycellH60};
%         % ====================== Inércia a Flexão ======================= %
%         L_w(3).D_I = {[7000 7000 7000]};
%         % ============= Fator Beta (define inércia a torsão) ============ %
%         L_w(3).beta = {[1 1 0.5]};
%         % ====================== Largura da Seção ======================= %
%         L_w(3).b = {[0.08 0.08 0.05]};


%% INICIALIZAÇÃO TIRANTE ANTERIOR %%
% ============================ INICIALIZAÇÃO ============================ %
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
%     %======================== Parâmetros de MEF ======================== %
%     L_w(7).cc = {};
%     L_w(7).idx_F = false;
%     L_w(7).connect = { 7 , 1 , @(L_w) L_w.posLong(1) , 'nearest' ; 7 , 2 , @(L_w) L_w.posLong(end) , 'nearest' };
%  
%     %======================== Parâmetros Tirante ======================= %
%         %======================= Fração de Corda ======================= %
%         L_w(7).fc = [0.3 0.3];
%        % ==================== Fração de Envergadura ==================== %
%         L_w(7).fb = [0.8 0.8];
%         
%    % ====================== Parâmetros Cicular ======================= %
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

    % ======================== Parâmetros de MEF ======================== %
    %L_w(5).cc = { 5 , 1 , 1:6 }';
    %L_w(5).idx_F = true;
%     L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
    %L_w(5).connect = {};

    % ======================== Parâmetros Globais ======================= %
        % ==================== Fração de Envergadura ==================== %
        %L_w(5).fb = [0.2 0.6 1];
        %L_w(5).fb0 = 0;
        % ======================= Fração de Corda ======================= %
        %L_w(5).fc =[0.7 0.7 0.7 0.7];
        %L_w(5).S_fc = [0.5 0.5 0.5];
        % ======================= Fração de Altura ====================== %
        %L_w(5).fh = [0 0 0 0];
        % ======================== Enflechamento ======================== %
        %L_w(5).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        %L_w(5).d_long = [nan nan nan];
        
    % ===================== Parâmetros Sanduiche I ====================== %
        % ====================== Altura do Núcleo ======================= %
        %L_w(5).H = [1 1 1 1];
        % ====================== Número de Camadas ====================== %
        %L_w(5).n_lam = [1 1;1 1; 1 1];
        % ====================== Fração com Camadas ===================== %
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

    % ======================== Parâmetros de MEF ======================== %
     %L_w(6).cc = { 6 , 1 , 2:2:6 }';
     %L_w(6).idx_F = true;
     %L_w(6).connect = { 6 , 5 , @(L_w) L_w.posLong(2) , 'nearest';6 , 5 , @(L_w) L_w.posLong(end) , 'nearest'};

    % ======================== Parâmetros Globais ======================= %
        % ==================== Fração de Envergadura ==================== %
         %L_w(6).fb = [0.4912 1]; %COLOCAR O FB DA MALHA DA ASA (Não é parametro da longarina)
         %L_w(6).fb0 = 0;
        % ======================= Fração de Corda ======================= %
         %L_w(6).fc = [0.1 0.1 0.1];
         %L_w(6).S_fc = [ 0.2 , 0.2];
        % ======================= Fração de Altura ====================== %
         %L_w(6).fh = [0 0 0];
        % =======================3 , 1 , @(L_w) L_w.posLong(5) , 'nearest'= Enflechamento ======================== %
         %L_w(6).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
         %L_w(6).d_long = [nan nan nan];
        
    % ======================== Parâmetros BA ============================ %
       % ======================= Fração de Corda ======================= %
         %L_w(6).fBA = [0.2 0.2 0.2];
        % ================ Número de Discretização em X ================= %
         %L_w(6).nx = 1000;
        % ====================== Número de Voltas ======================= %
         %L_w(6).n_voltas = [6 6;6 6];
        % ====================== Fração com Voltas ====================== %
         %L_w(6).f_voltas = [0.5 0.5];
% %% INICIALIZAÇÃO TUBINHO POSTERIOR %%
% % ============================ INICIALIZAÇÃO ============================ %
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
%     % ======================== Parâmetros de MEF ======================== %
%      L_w(5).cc = {5 , 1, 1:6}';
%      L_w(5).idx_F = true;
%      L_w(5).connect = {};
% 
%     % ======================== Parâmetros Globais ======================= %
%         % ==================== Fração de Envergadura ==================== %
%          L_w(5).fb = [1];
%          L_w(5).fb0 = 0;
%         % ======================= Fração de Corda ======================= %
%          L_w(5).fc = [ 0.11 0.11 ];
%          L_w(5).S_fc = [ nan ];
%         % ======================= Fração de Altura ====================== %
%          L_w(5).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%          L_w(5).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%          L_w(5).d_long = [nan nan nan];
%          % ======================= Parâmetros Tubinho ======================== %
%        % ====================== Diâmetro Interno ======================= %
%          L_w(5).D_interno = [2] * 1e-3;
%         % ====================== Diâmetro Externo ======================= %
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
%     % ======================== Parâmetros de MEF ======================== %
%      L_w(6).cc = {6 , 1, 1:6}';
%      L_w(6).idx_F = true;
%      L_w(6).connect = { 6 , 5 , @(L_w) L_w.posLong(end) , 'nearest'};
% 
%     % ======================== Parâmetros Globais ======================= %
%         % ==================== Fração de Envergadura ==================== %
%          L_w(6).fb = [1];
%          L_w(6).fb0 = 0;
%         % ======================= Fração de Corda ======================= %
%          L_w(6).fc = [ 0.31 0.31 ];
%          L_w(6).S_fc = [ nan ];
%         % ======================= Fração de Altura ====================== %
%          L_w(6).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%          L_w(6).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%          L_w(6).d_long = [nan nan nan];
%    % ======================= Parâmetros Tubinho ======================== %
%        % ====================== Diâmetro Interno ======================= %
%          L_w(6).D_interno = [2] * 1e-3;
%         % ====================== Diâmetro Externo ======================= %
%          L_w(6).D_externo = [3] * 1e-3;

% %% INICIALIZAÇÃO TUBINHO POSTERIOR %%
% % ============================ INICIALIZAÇÃO ============================ %
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
%     % ======================== Parâmetros de MEF ======================== %
%     L_w(7).cc = {}';
%     L_w(7).idx_F = true;
%     L_w(7).connect = { 7 , 3 , @(L_w) L_w.posLong(1) , 'nearest' };
% 
%     % ======================== Parâmetros Globais ======================= %
%         % ==================== Fração de Envergadura ==================== %
%         L_w(7).fb = [1];
%         L_w(7).fb0 = 1 - prod(1 - L_w(3).fb);
%         % ======================= Fração de Corda ======================= %
%         L_w(7).fc = [ 0.3 0.32 ];
%         L_w(7).S_fc = [ nan ];
%         % ======================= Fração de Altura ====================== %
%         L_w(7).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%         L_w(7).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%         L_w(7).d_long = [nan nan nan];
%         
%     % ======================= Parâmetros Tubinho ======================== %
%         % ====================== Diâmetro Interno ======================= %
%         L_w(7).D_interno = [6] * 1e-3;
%         % ====================== Diâmetro Externo ======================= %
%         L_w(7).D_externo = [8] * 1e-3;

end

