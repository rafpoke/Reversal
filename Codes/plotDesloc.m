function [] = plotDesloc( elem , nos , u , F_ext )
%% Carregando Valores %%
pontos = nos.pontos;
conectividade = elem.conect;

cc = nos.cc;

num_nos = size(pontos,1);

idx = 1:num_nos;

%% Corpo da Função %%
n_dim = 6;

% dx = u(([1:num_nos]-1)*n_dim+1);
% dy = u(([1:num_nos]-1)*n_dim+2);
% dz = u(([1:num_nos]-1)*n_dim+3);
% tx = u(([1:num_nos]-1)*n_dim+4);
% ty = u(([1:num_nos]-1)*n_dim+5);
% tz = u(([1:num_nos]-1)*n_dim+6);
% idx = 61:70;

dx = u((idx-1)*n_dim+1);
dy = u((idx-1)*n_dim+2);
dz = u((idx-1)*n_dim+3);% dz = dz * 0;
tx = u((idx-1)*n_dim+4);
ty = u((idx-1)*n_dim+5);
tz = u((idx-1)*n_dim+6);

ccx = cc((idx-1)*n_dim+1);
ccy = cc((idx-1)*n_dim+2);
ccz = cc((idx-1)*n_dim+3);
cctx = cc((idx-1)*n_dim+4);
ccty = cc((idx-1)*n_dim+5);
cctz = cc((idx-1)*n_dim+6);

%% Plot %%
% figure(1)
% clf
% C = dx;
% scatter3(nos(idx,1)+dx,nos(idx,2)+dy,nos(idx,3)+dz,'cdata',C)
% c = colorbar;
% c.Label.String = 'dX';
% hold on
% axis equal
% grid on
% 
% figure(2)
% clf
% C = dy;
% scatter3(pontos(idx,1)+dx,pontos(idx,2)+dy,pontos(idx,3)+dz,'cdata',C)
% c = colorbar;
% c.Label.String = 'Deslocamento em Y [m]';
% c.Label.FontSize = 14;
% hold on
% axis equal
% grid on
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% 
figure(3)
clf
C = dz;
scatter3(pontos(idx,1)+dx,pontos(idx,2)+dy,pontos(idx,3)+dz,'cdata',C)
c = colorbar;
c.Label.String = 'Deslocamento em Z [m]';
c.Label.FontSize = 14;
hold on
axis equal
grid on
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

% figure(4)
% clf
% C = tx;
% scatter3(nos(idx,1)+dx,nos(idx,2)+dy,nos(idx,3)+dz,'cdata',C)
% c = colorbar;
% c.Label.String = 'tX';
% hold on
% axis equal
% grid on
% 
% figure(5)
% clf
% C = ty;
% scatter3(nos(idx,1)+dx,nos(idx,2)+dy,nos(idx,3)+dz,'cdata',C)
% c = colorbar;
% c.Label.String = 'tY';
% hold on
% axis equal
% grid on
% 
% figure(6)
% clf
% C = tz;
% scatter3(nos(idx,1)+dx,nos(idx,2)+dy,nos(idx,3)+dz,'cdata',C)
% c = colorbar;
% c.Label.String = 'tZ';
% hold on
% axis equal
% grid on

% figure(7)
% clf
% C = double( ccz );
% scatter3(pontos(idx,1)+dx,pontos(idx,2)+dy,pontos(idx,3)+dz,'cdata',C)
% c = colorbar;
% c.Label.String = 'Engaste z';
% hold on
% axis equal
% grid on

for k = 1:size(conectividade,1)
%     i = idx(k);
    i = k;
    inicio = pontos(conectividade(i,1),:)+[dx(conectividade(i,1)),dy(conectividade(i,1)),dz(conectividade(i,1))];
    fim = pontos(conectividade(i,2),:)+[dx(conectividade(i,2)),dy(conectividade(i,2)),dz(conectividade(i,2))];
    linha = [inicio;fim];
    plot3(linha(:,1),linha(:,2),linha(:,3),'LineWidth',2,'Color',[0.4 0.4 0.4])
    hold on
    axis equal
    grid on
end

% h = scatter3(pontos(engaste(1),1)+dx(engaste(1)),pontos(engaste(1),2)+dy(engaste(1)),pontos(engaste(1),3)+dz(engaste(1)),'*','MarkerEdgeColor','r','LineWidth',4);
% legend(h,'Ponto de Engaste')

% plot3(nos(end-1,1)+dx(end-1),nos(end-1,2)+dy(end-1),nos(end-1,3)+dz(end-1),'*')
% for i = 1:num_nos
%     fx = F_ext((i-1)*n_dim+1);
%     fy = F_ext((i-1)*n_dim+2);
%     fz = F_ext((i-1)*n_dim+3);
%     
%     k = norm([fx,fy,fz]);
% % k = 1;
%     
%     quiver3(pontos(i,1)+dx(i),pontos(i,2)+dy(i),pontos(i,3)+dz(i),fx/k,fy/k,fz/k)
% % end
% 
% end

end

