function [ L ] = ETT_F_FEM( L , nos , F )
%% Carregando Valores %%
pontos = nos.pontos(:,2);
idx_long = nos.idx;
num_long = length(L);

if ~ isfield( L(1) , 'fem_F' )
    L(1).fem_F = [];
end


%% Corpo da Função %%
n_dim = 6;

for i = 1:num_long
    is_long = idx_long == i;
    idx_no1 = find( is_long == 1 , 1 , 'first' );
    idx_no2 = find( is_long == 1 , 1 , 'last' );
    idx_gl = ( (idx_no1-1)*n_dim + 1):( idx_no2*n_dim );

    L(i).fem_pontos = pontos( is_long );
    num_pontos = length( L(i).fem_pontos );
    
    if isempty( L(i).fem_F )
         L(i).fem_F = reshape( F( idx_gl ) , n_dim , num_pontos )';
    else
        L(i).fem_F(:,:,end+1) = reshape( F( idx_gl ) , n_dim , num_pontos )';
    end
    
end


end

