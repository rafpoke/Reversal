function [ geom_Asa ] = ETT_criaEH(aviao, ps, geom_Asa )
%% Inicializa as variáveis
c = aviao.spline_c_inf;
corda_raiz = ppval(c,0);
i_wl = aviao.i_w_inf;
dx_BF_BA_EH = aviao.dx_BF_BA_EH;

t = [1 aviao.t_EH];
offset_X = dx_BF_BA_EH + cosd(i_wl) * corda_raiz ;
l = [0 0];
perfis = repmat(aviao.perfil_EH, 1, 3);
i_wl_EH = aviao.i_w_EH;
twist = [0 0];

geom_Asa(3).x_BA = offset_X;                                                % [m] posição em x do BA da superfície (na raíz)
geom_Asa(3).h_BA = aviao.h_EH;                                              % [m] posição em z do BA da superfície (na raíz)
geom_Asa(3).c_sec = [1 t].*aviao.c_raiz_EH;                                 % [m] vetor de cordas da superfície
geom_Asa(3).b = aviao.b_EH;                                                 % [m] envergadura total da superfície
geom_Asa(3).l_BA = l;                                                       % [º] vetor de enflechamento no BA da superfície
geom_Asa(3).m_b = [aviao.fb_EH*aviao.b_EH/2 (1-aviao.fb_EH)*aviao.b_EH/2];  % [m] vetor de semi-envergaduras da superfície                         
geom_Asa(3).perfil = [22 22 23];                                                % [-] vetor de perfis da superfície
geom_Asa(3).i_w = i_wl_EH;                                                  % [º] incidência da superfície
geom_Asa(3).tw = twist;                                                     % [º] vetor de twists da superfície
% geom_Asa.sym = sym_w;                                                     % [bool] booleano de simetria da superfície                                              
% geom_Asa.rot = rot_w;                                                     % [bool] booleano de rotação da superfície (ex. empenagem vertical)       
geom_Asa(3).Malha = ETT_criaMalha( geom_Asa(3) , ps );                      % struct com a malha geométrica da superfície
% [ geom_Asa(3).Malha,~ ] = AeroMalhaEHGenetico( geom_eh, ps);
end                 