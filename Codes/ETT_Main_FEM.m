%% Função ETT_Main_FEM() %%
% ============================== Descrição ============================== %
% Essa função 


function [ Longarina ] = ETT_Main_FEM( Longarina , Cargas )

% for num_manobra = [1]
for num_manobra = 1:length(Cargas)
    F = arrayfun(@(x) x{1}', Cargas(num_manobra).F, 'UniformOutput', false);
    P = arrayfun(@(x) x{1}', Cargas(num_manobra).P, 'UniformOutput', false);
    F = [F{:}]';
    P = [P{:}]';
    
%     F = Cargas(num_manobra).F{2};%F = [0 0 100];
%     P = Cargas(num_manobra).P{2};%P = [0 1.8 0];

%     disp('Cuidado com ETT_Main_FEM!!! Tá só com a força da asa de cima!!!')
    % ====================== Programa de FEM principal ====================== %
    [ F_int , M_int , u , nos , elem , Longarina ] = ETT_LongarinaMEF( Longarina , F , P );
    % ==================== Atualização do vetor longarina =================== %
    [ Longarina ] = ETT_AtualizaLongarina( Longarina , nos , elem , F_int , M_int , u , num_manobra );
end
% ===================== Tensões sobre as longarinas ===================== %
Longarina = ETT_StressLongarina( Longarina );                                       % tensões atuantes sobre a longarina (Sy, Tyx, Tyz)
% ========================= Fator de segurança ========================== %
Longarina = ETT_FatorSeguranca( Longarina );  % fatores de segurança de cada região da longarina

end
