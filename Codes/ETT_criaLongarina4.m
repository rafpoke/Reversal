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


function [ L_w ] = ETT_criaLongarina4( materiais )
%% INICIALIZAÇÃO LONGARINA ASA %%
% ============================ INICIALIZAÇÃO ============================ %
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

    % ======================== Parâmetros de MEF ======================== %
    L_w(1).cc = { 1 , 1 , 1:6 }';
    L_w(1).idx_F = true;
    L_w(1).connect = {};

    % ======================== Parâmetros Globais ======================= %
        % ==================== Fração de Envergadura ==================== %
        L_w(1).fb = [0.5 1];
        L_w(1).fb0 = 0;
        % ======================= Fração de Corda ======================= %
        L_w(1).fc = [0.4 0.4 0.2];
        L_w(1).S_fc = [ 0.5 , 0.2 ];
        % ======================= Fração de Altura ====================== %
        L_w(1).fh = [0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(1).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(1).d_long = [nan nan nan];
        
    % ======================== Parâmetros BA ============================ %
        % ======================= Fração de Corda ======================= %
        L_w(1).fBA = [0.2 0.1 0.05];
        % ================ Número de Discretização em X ================= %
        L_w(1).nx = 100;
        % ====================== Número de Voltas ======================= %
        L_w(1).n_voltas = [1 1;1 1];
        % ====================== Fração com Voltas ====================== %
        L_w(1).f_voltas = [0.5 0.5];

%% INICIALIZAÇÃO LONGARINA AUXILIAR %%
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

    % ======================== Parâmetros de MEF ======================== %
    L_w(2).cc = { 2 , 1 , 1:6 }';
    L_w(2).idx_F = true;
    L_w(2).connect = { 2 , 1 , @(L_w) L_w.posLong(end) , 'nearest' };
%     L_w(2).connect = {};

    % ======================== Parâmetros Globais ======================= %
        % ==================== Fração de Envergadura ==================== %
        L_w(2).fb = [0.22];
        L_w(2).fb0 = 0;
        % ======================= Fração de Corda ======================= %
        L_w(2).fc = [0.2 0.2];
        L_w(2).S_fc = 0.5;
        % ======================= Fração de Altura ====================== %
        L_w(2).fh = [0 0];
        % ======================== Enflechamento ======================== %
        L_w(2).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(2).d_long = [nan nan nan];
        
    % ======================= Parâmetros Circular ======================= %
        % ========================== Diâmetro =========================== %
        L_w(2).D = [9/8] * 25.4e-3;
        % ====================== Número de Voltas ======================= %
        L_w(2).n_voltas = [3 2];
        % ====================== Fração com Voltas ====================== %
        L_w(2).f_voltas = [0.5];
        
%% INICIALIZAÇÃO LONGARINAS - ASA SUPERIOR %%
% ============================ INICIALIZAÇÃO ============================ %
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

    % ======================== Parâmetros de MEF ======================== %
    L_w(3).cc = { 1 , 1 , 1:6 }';
    L_w(3).idx_F = true;
    L_w(3).connect = {};

    % ======================== Parâmetros Globais ======================= %
        % ==================== Fração de Envergadura ==================== %
        L_w(3).fb = [0.5 1];
        L_w(3).fb0 = 0;
        % ======================= Fração de Corda ======================= %
        L_w(3).fc = [0.4 0.4 0.2];
        L_w(3).S_fc = [ 0.5 , 0.2 ];
        % ======================= Fração de Altura ====================== %
        L_w(3).fh = [0 0 0];
        % ======================== Enflechamento ======================== %
        L_w(3).l_long = [nan nan nan];
        % =========================== Diedro ============================ %
        L_w(3).d_long = [nan nan nan];
        
    % ======================== Parâmetros BA ============================ %
        % ======================= Fração de Corda ======================= %
        L_w(3).fBA = [0.3 0.1 0.05];
        % ================ Número de Discretização em X ================= %
        L_w(3).nx = 100;
        % ====================== Número de Voltas ======================= %
        L_w(3).n_voltas = [1 1;1 1];
        % ====================== Fração com Voltas ====================== %
        L_w(3).f_voltas = [0.5 0.5];

        
%% INICIALIZAÇÃO TIRANTE POSTERIOR %%
% ============================ INICIALIZAÇÃO ============================ %
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

    % ======================== Parâmetros de MEF ======================== %
   % L_w(4).cc = {};
   % L_w(4).idx_F = false;
   % L_w(4).connect = { 4 , 2 , @(L_w) L_w.posLong(1) , 'nearest' ; 4 , 3 , @(L_w) L_w.posLong(end) , 'nearest' };

    % ======================== Parâmetros Tirante ======================= %
        % ======================= Fração de Corda ======================= %
       % L_w(4).fc = [0.15 0.15];
        % ==================== Fração de Envergadura ==================== %
       % L_w(4).fb = [0 0];
        % ======================= Número de seções ====================== %
       % L_w(4).num_T = 5;
    % ======================== Parâmetros Treliça ======================= %
        % =========================== Largura =========================== %
       % L_w(4).b = [0.08 0.08;0.08 0.08];
        % ==================== Curvatura da Largura ===================== %
       % L_w(4).Sb = [0.8642 0.8642];
        % ======================= Comprimento em X ====================== %
       % L_w(4).h = [0.08 0.04 0.08];
        % =================== Curvatura do Comprimento ================== %
       % L_w(4).Sh = 0.92;


%% INICIALIZAÇÃO TIRANTE ANTERIOR %%
% ============================ INICIALIZAÇÃO ============================ %
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
%     % ======================== Parâmetros de MEF ======================== %
%     L_w(5).cc = {};
%     L_w(5).idx_F = false;
%     L_w(5).connect = { 5 , 1 , @(L_w) L_w.posLong(1) , 'nearest' ; 5 , 3 , @(L_w) L_w.posLong(end) , 'nearest' };
% 
%     % ======================== Parâmetros Tirante ======================= %
%         % ======================= Fração de Corda ======================= %
%         L_w(5).fc = [0.45 0.3];
%         % ==================== Fração de Envergadura ==================== %
%         L_w(5).fb = [0 0];
%         
%     % ====================== Parâmetros Sanduíche ======================= %
%         % ======================= Material Núcleo ======================= %
%         L_w(5).material_c = {materiais.divinycell};
%         % ====================== Inércia a Flexão ======================= %
%         L_w(5).D_I = {[600 600 0.5]};
%         % ============= Fator Beta (define inércia a torsão) ============ %
%         L_w(5).beta = {[1 1 0.5];[1 1 0.5]};
%         % ====================== Largura da Seção ======================= %
%         L_w(5).b = {[0.05 0.05 0.5]};
%         
% %% INICIALIZAÇÃO TUBINHO POSTERIOR %%
% % ============================ INICIALIZAÇÃO ============================ %
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
%     % ======================== Parâmetros de MEF ======================== %
%     L_w(6).cc = {}';
%     L_w(6).idx_F = true;
%     L_w(6).connect = { 6 , 1 , @(L_w) L_w.posLong(1) , 'nearest' };
% 
%     % ======================== Parâmetros Globais ======================= %
%         % ==================== Fração de Envergadura ==================== %
%         L_w(6).fb = [1];
%         L_w(6).fb0 = 1 - prod(1 - L_w(1).fb);
%         % ======================= Fração de Corda ======================= %
%         L_w(6).fc = [ 0.355 0.32 ];
%         L_w(6).S_fc = [ nan ];
%         % ======================= Fração de Altura ====================== %
%         L_w(6).fh = [0 0];
%         % ======================== Enflechamento ======================== %
%         L_w(6).l_long = [nan nan nan];
%         % =========================== Diedro ============================ %
%         L_w(6).d_long = [nan nan nan];
%         
%     % ======================= Parâmetros Tubinho ======================== %
%         % ====================== Diâmetro Interno ======================= %
%         L_w(6).D_interno = [2] * 1e-3;
%         % ====================== Diâmetro Externo ======================= %
%         L_w(6).D_externo = [4] * 1e-3;
% 
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

