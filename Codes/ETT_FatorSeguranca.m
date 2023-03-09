function [ Longarina ] = ETT_FatorSeguranca( Longarina )
%% CARREGANDO VALORES %%
% Longarina = aviao.Longarina;
material = {Longarina.material};
Sx = {Longarina.Sx};
T = {Longarina.T};
posLong = {Longarina.posLong};

%% CORPO DA FUNÇÃO %%
for i = 1:length(Longarina)-2
    Longarina(i).CS = nan(1,size(Sx{i},2));
    Longarina(i).CS_c = nan(1,size(Sx{i},2));
    for j = 1:length(Longarina(i).config)
        idx = posLong{i}(j):posLong{i}(j+1) - 1;
        
        if strcmp(material{i}{j}.criterio,'Tsai-Wu') && (strcmp(Longarina(i).config{j},'sanduiche')||strcmp(Longarina(i).config{j},'sanduiche_I')||strcmp(Longarina(i).config{j},'sanduiche_H'))        
            Sx_c = {Longarina.Sx_c};
            T_c = {Longarina.T_c};
            material_c = {Longarina.material_c};
                   
            Longarina(i).CS(idx) = ETT_FatorSegurancaSimples(material{i}{j},Sx{i}(:,idx,:),T{i}(:,idx,:));
            Longarina(i).CS_c(idx) = ETT_FatorSegurancaSimples(material_c{i}{j},Sx_c{i}(:,idx,:),T_c{i}(:,idx,:));
            
        elseif strcmp(material{i}{j}.criterio,'Tsai-Wu')
            Longarina(i).CS(idx) = ETT_FatorSegurancaSimples(material{i}{j},Sx{i}(:,idx,:),T{i}(:,idx,:));            
        end     
    end
end

end

%% Função ETT_FatorSegurancaSimples() %%
function CS = ETT_FatorSegurancaSimples(material,Sx,T)
%% CRITÉRIO DE TSAI-WU (ESTADO PLANO DE TENSÃO) %%
% ======================= CARREGANDO VALORES ======================== %
F1t = material.St1;
F1c = material.Sc1;
F2t = material.St2;
F2c = material.Sc2;
F6 = material.T1;

% ===================== COEFICIENTES DO CRITÉRIO ==================== %
f1 = 1/F1t - 1/F1c;
f11 = 1/(F1t*F1c);
f2 = 1/F2t - 1/F2c;
f22 = 1/(F2t*F2c);
f66 = 1/F6^2;
f12 = -(1/2)*sqrt(f11*f22);
% ======================== FATOR DE SEGURANÇA ======================= %
C1 = f1*Sx;
C2 = f11*Sx.^2 + f66*T.^2;
R = sqrt( abs( C1.^2 + 4*C2 ) );
CS = (R - C1)./(2*C2);
CS = min( min(CS,[],1) ,[],3 );

% if strcmp(material.criterio,'von Mises')
%     
% end
end
