function [] = tableElem(elem)

A = elem.A';
Ix = elem.Ix';
Iy = elem.Iy';
Iz = elem.Iz';
E = elem.E';
G = elem.G';
L = elem.L';

idx = elem.idx;
num = [1:length(A)]';

table(num,idx,A,Ix,Iy,Iz,E,G,L)

end

