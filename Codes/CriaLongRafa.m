Long_25_05=struct;
%%%%%%Longarina Inferior%%%%%%%%%%%
Long_25_05.L_inf.E_unidade = "Pa";
Long_25_05.L_inf.G_unidade = "Pa";
Long_25_05.L_inf.Inercias_unidade = "m^4";
Long_25_05.L_inf.E_nucleo = L_w(1).En;
Long_25_05.L_inf.G_nucleo = L_w(1).Gn;
Long_25_05.L_inf.E_face = L_w(1).Ex;
Long_25_05.L_inf.G_face = L_w(1).Gx;
Long_25_05.L_inf.Ix = L_w(1).Ix;
Long_25_05.L_inf.Iy = L_w(1).Iy;
Long_25_05.L_inf.IZ = L_w(1).Iz;
Long_25_05.L_inf.EIy = L_w(1).EIy;
Long_25_05.L_inf.EIz = L_w(1).EIz;
Long_25_05.L_inf.GJx = L_w(1).GJx;
%%%%%%% Longarina Superior %%%%%%%%%5
Long_25_05.L_sup.E_unidade = "Pa";
Long_25_05.L_sup.G_unidade = "Pa";
Long_25_05.L_sup.Inercias_unidade = "m^4";
Long_25_05.L_sup.E_nucleo = L_w(2).En;
Long_25_05.L_sup.G_nucleo = L_w(2).Gn;
Long_25_05.L_sup.E_face = L_w(2).Ex;
Long_25_05.L_sup.G_face = L_w(2).Gx;
Long_25_05.L_sup.Ix = L_w(2).Ix;
Long_25_05.L_sup.Iy = L_w(2).Iy;
Long_25_05.L_sup.IZ = L_w(2).Iz;
Long_25_05.L_sup.EIy = L_w(2).EIy;
Long_25_05.L_sup.EIz = L_w(2).EIz;
Long_25_05.L_sup.GJx = L_w(2).GJx;

save("Long_25_05.mat","Long_25_05")