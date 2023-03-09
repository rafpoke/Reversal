%% Fun��o ETT_criaGeom() %%
% ============================== Descri��o ============================== %
% Essa fun��o cria struct que descrevem a geometria de cada superf�cie
% aerodin�mica. Empregando essa abordagem � poss�vel generalizar as fun��es
% respons�veis por calcular a massa desas superf�cies.

function [ geom_Asa ] = ETT_criaGeom2( aviao , ps )
%% CARREGANDO VALORES %%
% ================================= Asa ================================= %
% h_BA = aviao.h_BA;
% h_BA = aviao.h_asa1_BA;
% c_sec = aviao.c_sec;
% b_total = aviao.b;
% l_BA = aviao.l_BA;
% m_b = aviao.m_b;
% perfil = aviao.perfil;
% % i_w = aviao.i_w;
% i_w = aviao.i_w1;
% tw_w = aviao.tw;
% x_BA_w = 0;
% sym_w = true;
% ================================== EH ================================= %
% h_BA_h = aviao.h_BA_h;
% c_sec_h = aviao.c_sec_h;
% b_total_h = aviao.b_h;
% % l_BA_h = aviao.l_BA_h;
% l_BA_h = 0;
% m_b_h = b_total_h/2;                                                        % Assumo que a EH tem apenas 1 se��o.
% perfil_h = [aviao.perfil_h aviao.perfil_h];                                 % Assumo que EH tem apenas 1 perfil
% i_h = aviao.i_h;
% tw_h = 0;
% x_BA_h = aviao.d_BA_x;
% sym_h = true;
% ================================== EV ================================= %
% h_BA_v = 0;
% c_sec_v = [0.3 0.15];
% b_total_v = 0.8;
% l_BA_v = 0;
% m_b_v = b_total_v/2;                                                          % Assumo que a EH tem apenas 1 se��o.
% perfil_v = [7 7];                                                           % Assumo que EH tem apenas 1 perfil
% i_v = 0;                                                                    % Assumo que a EV n�o tem incid�ncia
% tw_v = 0;
% x_BA_v = 0;
% sym_v = false;

%% CRIANDO STRUCTS %%
% =============================== geom_Asa ============================== %
geom_Asa.xBA = 0;                                                           % [m] posi��o em x do BA da superf�cie (na ra�z)
geom_Asa.zBA = aviao.h_BA_inf;                                              % [m] posi��o em z do BA da superf�cie (na ra�z)
geom_Asa.spline_c = aviao.spline_c_inf;                                     % [m] vetor de cordas da superf�cie
geom_Asa.b = aviao.b;                                                       % [m] envergadura total da superf�cie
geom_Asa.b_sec = aviao.b_sec_inf;
geom_Asa.spline_l = aviao.spline_l_inf;                                     % [�] vetor de enflechamento no BA da superf�cie
geom_Asa.spline_tw = aviao.spline_tw_inf;                                   % [m] vetor de semi-envergaduras da superf�cie                         
geom_Asa.perfil = aviao.perfil_1;                                             % [-] vetor de perfis da superf�cie
geom_Asa.i_w = aviao.i_w_inf;                                               % [�] incid�ncia da superf�cie
% geom_Asa.sym = sym_w;                                                       % [bool] booleano de simetria da superf�cie                                              
% geom_Asa.rot = rot_w;                                                       % [bool] booleano de rota��o da superf�cie (ex. empenagem vertical)       
geom_Asa.Malha = ETT_criaMalha2( geom_Asa , ps );                            % struct com a malha geom�trica da superf�cie

% =============================== geom_Asa ============================== %
geom_Asa(2).xBA = 0;                                                           % [m] posi��o em x do BA da superf�cie (na ra�z)
geom_Asa(2).zBA = aviao.h_BA_sup;                                              % [m] posi��o em z do BA da superf�cie (na ra�z)
geom_Asa(2).spline_c = aviao.spline_c_sup;                                     % [m] vetor de cordas da superf�cie
geom_Asa(2).b = aviao.b;                                                       % [m] envergadura total da superf�cie
geom_Asa(2).b_sec = aviao.b_sec_sup;
geom_Asa(2).spline_l = aviao.spline_l_sup;                                     % [�] vetor de enflechamento no BA da superf�cie
geom_Asa(2).spline_tw = aviao.spline_tw_sup;                                   % [m] vetor de semi-envergaduras da superf�cie                         
geom_Asa(2).perfil = aviao.perfil_2;                                             % [-] vetor de perfis da superf�cie
geom_Asa(2).i_w = aviao.i_w_sup;                                               % [�] incid�ncia da superf�cie
% geom_Asa.sym = sym_w;                                                       % [bool] booleano de simetria da superf�cie                                              
% geom_Asa.rot = rot_w;                                                       % [bool] booleano de rota��o da superf�cie (ex. empenagem vertical)       
geom_Asa(2).Malha = ETT_criaMalha2( geom_Asa(2) , ps );                            % struct com a malha geom�trica da superf�cie

end

