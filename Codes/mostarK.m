function [] = mostarK(K_loc,K_glob,struct_fem)
conectividade = struct_fem.conectividade;
n_dim = 6;
num_elem = size(conectividade,1);

for i = 1:num_elem
    no1 = conectividade(i,1);
    no2 = conectividade(i,2);
    idx1 = [(no1-1)*n_dim+1:no1*n_dim];
    idx2 = [(no2-1)*n_dim+1:no2*n_dim];
    
    str = 'Elemento ' + string(i);
    disp(str)
    disp('Matriz k' + string(no1) + string(no1))
    disp(K_loc(1:n_dim,1:n_dim,i));
    disp('Matriz k' + string(no2) + string(no2))
    disp(K_loc(n_dim+1:2*n_dim,n_dim+1:2*n_dim,i));
    disp('Matriz k' + string(no1) + string(no2))
    disp(K_loc(1:n_dim,n_dim+1:2*n_dim,i));
end
end

