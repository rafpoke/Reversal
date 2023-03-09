%% Função ETT_AtualizaLongarina() %%

function [ L ] = ETT_AtualizaLongarina( Longarina , nos , elem , F_int , M_int , u , num_manobra )
%% Carregando Valores %%
L = Longarina;
idx_long_nos = nos.idx;
idx_long_elem = elem.idx;
num_long = length(L);
conect = elem.conect;
nos_pontos = nos.pontos;

%% Atualização da Longarina %%
% ===================== Ordenação dos Deslocamentos ===================== %
n_dim = 6;
u = reshape(u,n_dim,[]);

% ========================= Pontos dos Elementos ======================== %
num_elem = size(conect,1);
aux = nos_pontos(conect',:);
elem_pontos = (aux((1:num_elem)*2,:)+aux((1:num_elem)*2-1,:))/2;

% ============================= Atualização ============================= %
for i = 1:num_long
    L(i).F(:,:,num_manobra) = F_int( idx_long_elem == i , : );
    L(i).M(:,:,num_manobra) = M_int( idx_long_elem == i , : );    
    L(i).u(:,:,num_manobra) = u( : , idx_long_nos == i );
    
    L(i).elem_pontos = elem_pontos( idx_long_elem == i , : );
end

end

