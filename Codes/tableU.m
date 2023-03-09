function [] = tableU(u,F_int,M_int,F_ext,elem_struct)
n_dim = 6;
num_nos = length(u)/n_dim;
num_elem = size(F_int,1);

idx = elem_struct.idx;

no = cell(num_nos,1);
elem = cell(num_elem,1);
x = u(([1:num_nos]-1)*n_dim+1);
y = u(([1:num_nos]-1)*n_dim+2);
z = u(([1:num_nos]-1)*n_dim+3);
tx = u(([1:num_nos]-1)*n_dim+4);
ty = u(([1:num_nos]-1)*n_dim+5);
tz = u(([1:num_nos]-1)*n_dim+6);

Fx = F_int(:,1);
Fy = F_int(:,2);
Fz = F_int(:,3);
Mx = M_int(:,1);
My = M_int(:,2);
Mz = M_int(:,3);

F1 = F_ext(([1:num_nos]-1)*n_dim+1);
F2 = F_ext(([1:num_nos]-1)*n_dim+2);
F3 = F_ext(([1:num_nos]-1)*n_dim+3);
M1 = F_ext(([1:num_nos]-1)*n_dim+4);
M2 = F_ext(([1:num_nos]-1)*n_dim+5);
M3 = F_ext(([1:num_nos]-1)*n_dim+6);


for i = 1:num_nos
    no{i} = string(i);
end

for i = 1:num_elem
    elem{i} = string(i);
end

table(no,x,y,z,tx,ty,tz)
table(elem,idx,Fx,Fy,Fz,Mx,My,Mz)
table(no,F1,F2,F3,M1,M2,M3)

a = [F1,F2,F3,M1,M2,M3];
assignin('base','a',a)
end

