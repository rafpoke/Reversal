%% Função ETT_criaMateriais() %%
function [materiais] = ETT_criaMateriais()
% ============================ INICIALIZAÇÃO ============================ %
materiais = struct;

%% MATERIAIS %%
% =============================== CARBONO =============================== %
carbono = struct;
carbono.rho = 1600;                                                         % densidade [kg/m^3]
carbono.E1 = 70e9;                                                          % módulo de elasticidade na direção principal [Pa]
carbono.E2 = 70e9;                                                          % módulo de elasticidade na direção ortogonal [Pa]
carbono.G1 = 19.2e9;                                                           % módulo de cisalhamento na direção principal [Pa]
carbono.G2 = 19.2e9;                                                           % módulo de cisalhamento na direção ortogonal [Pa]
carbono.St1 = 600e6;                                                        % tensão de falha por tração na direção principal [Pa]
carbono.St2 = 600e6;                                                        % tensão de falha por tração na direção ortogonal [Pa]
carbono.Sc1 = 570e6;                                                        % tensão de falha por compressão na direção principal [Pa]
carbono.Sc2 = 570e6;                                                        % tensão de falha por compressão na direção ortogonal [Pa]
carbono.T1 = 60e6;                                                          % tensão de falha por cisalhamento na direção principal [Pa]
carbono.T2 = 60e6;                                                          % tensão de falha por cisalhamento na direção ortogonal [Pa]
carbono.criterio = 'Tsai-Wu';                                               % critério de falha adotado para o material

materiais.carbono = carbono;
% =========================== CARBONO-ARAMIDA =========================== %
ck = struct;
ck.rho = 1130;                                                              % densidade [kg/m^3]
ck.E1 = 32.3e9;                                                             % módulo de elasticidade na direção principal [Pa]
ck.E2 = 32.3e9;                                                             % módulo de elasticidade na direção ortogonal [Pa]
ck.G1 = 12.8e9;                                                             % módulo de cisalhamento na direção principal [Pa]
ck.G2 = 12.8e9;                                                             % módulo de cisalhamento na direção ortogonal [Pa]
ck.St1 = 450e6;                                                             % tensão de falha por tração na direção principal [Pa]
ck.St2 = 450e6;                                                             % tensão de falha por tração na direção ortogonal [Pa]
ck.Sc1 = 400e6;                                                             % tensão de falha por compressão na direção principal [Pa]
ck.Sc2 = 400e6;                                                             % tensão de falha por compressão na direção ortogonal [Pa]
ck.T1 = 45e6;                                                               % tensão de falha por cisalhamento na direção principal [Pa]
ck.T2 = 45e6;                                                               % tensão de falha por cisalhamento na direção ortogonal [Pa]
ck.criterio = 'Tsai-Wu';                                                    % critério de falha adotado para o material

materiais.ck = ck;
% =============================== BALSA 1A ============================== %
balsa_1A = struct;
balsa_1A.rho = 299;                                                         % densidade [kg/m^3]
balsa_1A.E1 = 1989e6;                                                       % módulo de elasticidade na direção principal [Pa]
balsa_1A.E2 = 91e6;                                                         % módulo de elasticidade na direção ortogonal [Pa]
balsa_1A.G1 = 107e6;                                                        % módulo de cisalhamento na direção principal [Pa]
balsa_1A.G2 = 74e6;                                                         % módulo de cisalhamento na direção ortogonal [Pa]
balsa_1A.St1 = 44.1e6;                                                      % tensão de falha por tração na direção principal [Pa]
balsa_1A.St2 = 2.5e6;                                                       % tensão de falha por tração na direção ortogonal [Pa]
balsa_1A.Sc1 = 19.8e6;                                                      % tensão de falha por compressão na direção principal [Pa]
balsa_1A.Sc2 = 1.6e6;                                                       % tensão de falha por compressão na direção ortogonal [Pa]
balsa_1A.T1 = 5.3e6;                                                        % tensão de falha por cisalhamento na direção principal [Pa]
balsa_1A.T2 = 5.3e6;                                                        % tensão de falha por cisalhamento na direção ortogonal [Pa]
balsa_1A.criterio = 'Tsai-Wu';                                              % critério de falha adotado para o material

materiais.balsa_1A = balsa_1A;
% =============================== BALSA 2A ============================== %
balsa_2A = struct;
balsa_2A.rho = 299;                                                         % densidade [kg/m^3]
balsa_2A.E1 = 1230e6;                                                       % módulo de elasticidade na direção principal [Pa]
balsa_2A.E2 = 57e6;                                                         % módulo de elasticidade na direção ortogonal [Pa]
balsa_2A.G1 = 66e6;                                                         % módulo de cisalhamento na direção principal [Pa]
balsa_2A.G2 = 46e6;                                                         % módulo de cisalhamento na direção ortogonal [Pa]
balsa_2A.St1 = 33.3e6;                                                      % tensão de falha por tração na direção principal [Pa]
balsa_2A.St2 = 1.9e6;                                                       % tensão de falha por tração na direção ortogonal [Pa]
balsa_2A.Sc1 = 15.0e6;                                                      % tensão de falha por compressão na direção principal [Pa]
balsa_2A.Sc2 = 1.2e6;                                                       % tensão de falha por compressão na direção ortogonal [Pa]
balsa_2A.T1 = 4.0e6;                                                        % tensão de falha por cisalhamento na direção principal [Pa]
balsa_2A.T2 = 4.0e6;                                                        % tensão de falha por cisalhamento na direção ortogonal [Pa]
balsa_2A.criterio = 'Tsai-Wu';                                              % critério de falha adotado para o material

materiais.balsa_2A = balsa_2A;
% =============================== BALSA 3A ============================== %
balsa_3A = struct;
balsa_3A.rho = 207;                                                         % densidade [kg/m^3]
balsa_3A.E1 = 990e6;                                                        % módulo de elasticidade na direção principal [Pa]
balsa_3A.E2 = 46e6;                                                         % módulo de elasticidade na direção ortogonal [Pa]
balsa_3A.G1 = 53e6;                                                         % módulo de cisalhamento na direção principal [Pa]
balsa_3A.G2 = 37e6;                                                         % módulo de cisalhamento na direção ortogonal [Pa]
balsa_3A.St1 = 20.7e6;                                                      % tensão de falha por tração na direção principal [Pa]
balsa_3A.St2 = 1.2e6;                                                       % tensão de falha por tração na direção ortogonal [Pa]
balsa_3A.Sc1 = 9.3e6;                                                       % tensão de falha por compressão na direção principal [Pa]
balsa_3A.Sc2 = 0.7e6;                                                       % tensão de falha por compressão na direção ortogonal [Pa]
balsa_3A.T1 = 2.5e6;                                                        % tensão de falha por cisalhamento na direção principal [Pa]
balsa_3A.T2 = 2.5e6;                                                        % tensão de falha por cisalhamento na direção ortogonal [Pa]
balsa_3A.criterio = 'Tsai-Wu';                                              % critério de falha adotado para o material

materiais.balsa_3A = balsa_3A;
% =============================== BALSA 4A ============================== %
balsa_4A = struct;
balsa_4A.rho = 151;                                                         % densidade [kg/m^3]
balsa_4A.E1 = 638e6;                                                        % módulo de elasticidade na direção principal [Pa]
balsa_4A.E2 = 29e6;                                                         % módulo de elasticidade na direção ortogonal [Pa]
balsa_4A.G1 = 34e6;                                                         % módulo de cisalhamento na direção principal [Pa]
balsa_4A.G2 = 24e6;                                                         % módulo de cisalhamento na direção ortogonal [Pa]
balsa_4A.St1 = 7.3e6;                                                       % tensão de falha por tração na direção principal [Pa]
balsa_4A.St2 = 0.4e6;                                                       % tensão de falha por tração na direção ortogonal [Pa]
balsa_4A.Sc1 = 3.3e6;                                                       % tensão de falha por compressão na direção principal [Pa]
balsa_4A.Sc2 = 0.3e6;                                                       % tensão de falha por compressão na direção ortogonal [Pa]
balsa_4A.T1 = 0.9e6;                                                        % tensão de falha por cisalhamento na direção principal [Pa]
balsa_4A.T2 = 0.9e6;                                                        % tensão de falha por cisalhamento na direção ortogonal [Pa]
balsa_4A.criterio = 'Tsai-Wu';                                              % critério de falha adotado para o material

materiais.balsa_4A = balsa_4A;
% ============================= BALSA EXTRA ============================= %
balsa_X = struct;
balsa_X.rho = 91;                                                           % densidade [kg/m^3]
balsa_X.E1 = 374e6;                                                         % módulo de elasticidade na direção principal [Pa]
balsa_X.E2 = 17e6;                                                          % módulo de elasticidade na direção ortogonal [Pa]
balsa_X.G1 = 20e6;                                                          % módulo de cisalhamento na direção principal [Pa]
balsa_X.G2 = 14e6;                                                          % módulo de cisalhamento na direção ortogonal [Pa]
balsa_X.St1 = 5.8e6;                                                        % tensão de falha por tração na direção principal [Pa]
balsa_X.St2 = 0.3e6;                                                        % tensão de falha por tração na direção ortogonal [Pa]
balsa_X.Sc1 = 2.6e6;                                                        % tensão de falha por compressão na direção principal [Pa]
balsa_X.Sc2 = 0.2e6;                                                        % tensão de falha por compressão na direção ortogonal [Pa]
balsa_X.T1 = 0.7e6;                                                         % tensão de falha por cisalhamento na direção principal [Pa]
balsa_X.T2 = 0.7e6;                                                         % tensão de falha por cisalhamento na direção ortogonal [Pa]
balsa_X.criterio = 'Tsai-Wu';                                               % critério de falha adotado para o material

materiais.balsa_X = balsa_X;
% ========================= COMPENSADO DE TAUARÍ ======================== %
tauari = struct;
tauari.rho = 554;                                                           % densidade [kg/m^3]
tauari.E1 = 1273e6;                                                         % módulo de elasticidade na direção principal [Pa]
tauari.E2 = 1690e6;                                                         % módulo de elasticidade na direção ortogonal [Pa]
tauari.G1 = 116e6;                                                          % módulo de cisalhamento na direção principal [Pa]
tauari.G2 = 116e6;                                                          % módulo de cisalhamento na direção ortogonal [Pa]
tauari.St1 = 27.9e6;                                                        % tensão de falha por tração na direção principal [Pa]
tauari.St2 = 32.3e6;                                                        % tensão de falha por tração na direção ortogonal [Pa]
tauari.Sc1 = 40.3e6;                                                        % tensão de falha por compressão na direção principal [Pa]
tauari.Sc2 = 22.4e6;                                                        % tensão de falha por compressão na direção ortogonal [Pa]
tauari.T1 = 9.5e6;                                                          % tensão de falha por cisalhamento na direção principal [Pa]
tauari.T2 = 9.5e6;                                                          % tensão de falha por cisalhamento na direção ortogonal [Pa]
tauari.criterio = 'Tsai-Wu';                                                % critério de falha adotado para o material

materiais.tauari = tauari;
% ============================= DIVINYCELL ============================== %
divinycell = struct;
divinycell.rho = 60;
divinycell.E1 = 0.075e9;
divinycell.E2 = 0.075e9;
divinycell.G1 = 0.02e9;
divinycell.G2 = 0.02e9;
divinycell.St1 = 1.8e6;
divinycell.St2 = 1.8e6;
divinycell.Sc1 = 0.9e6;
divinycell.Sc2 = 0.9e6;
divinycell.T1 = 0.76e6;
divinycell.T2 = .76e6;
divinycell.criterio = 'Tsai-Wu';

materiais.divinycell = divinycell;
% ============================= ROHACELL_31A ============================== %
rohacell_31A = struct;
rohacell_31A.rho = 32;
rohacell_31A.E1 = 0.036e9;
rohacell_31A.E2 = 0.036e9;
rohacell_31A.G1 = 0.013e9;
rohacell_31A.G2 = 0.013e9;
rohacell_31A.St1 = 1.0e6;
rohacell_31A.St2 = 1.0e6;
rohacell_31A.Sc1 = 0.4e6;
rohacell_31A.Sc2 = 0.4e6;
rohacell_31A.T1 = 0.4e6;
rohacell_31A.T2 = 0.4e6;
rohacell_31A.criterio = 'Tsai-Wu';

materiais.rohacell_31A = rohacell_31A;
% ============================== ALUMÍNIO =============================== %
aluminio = struct;
aluminio.isotropico = false;                                                % material não isotrópico
aluminio.rho = 2810;                                                        % densidade [kg/m^3]
aluminio.E1 = 71.1e9;                                                       % módulo de elasticidade na direção principal [Pa]
aluminio.E2 = 71.18e9;                                                      % módulo de elasticidade na direção ortogonal [Pa]
aluminio.G1 = 26.9e9;                                                       % módulo de cisalhamento na direção principal [Pa]
aluminio.G2 = 26.9e9;                                                       % módulo de cisalhamento na direção ortogonal [Pa]
aluminio.St1 = 572e6;                                                       % tensão de falha por tração na direção principal [Pa]
aluminio.St2 = 572e6;                                                       % tensão de falha por tração na direção ortogonal [Pa]
aluminio.Sc1 = 572e6;                                                       % tensão de falha por compressão na direção principal [Pa]
aluminio.Sc2 = 572e6;                                                       % tensão de falha por compressão na direção ortogonal [Pa]
aluminio.T1 = 331e6;                                                        % tensão de falha por cisalhamento na direção principal [Pa]
aluminio.T2 = 331e6;                                                        % tensão de falha por cisalhamento na direção ortogonal [Pa]
aluminio.criterio = 'Von Mises';                                            % critério de falha adotado para o material

materiais.aluminio = aluminio;

end

