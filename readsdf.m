tic
clear all
% figure
i   = 0;
%E_m = 4e12;
% matlabpool local 4
for n =59
    if floor(log10(n)) == 0
    %deta_x = 80/2400*1e-6;
%     deta_y = 160/1600*1e-6;
   % load(sprintf('ele%d.mat',n))
        data=GetDataSDF(sprintf('000%d.sdf',n));
    elseif floor(log10(n)) == 1
        data=GetDataSDF(sprintf('00%d.sdf',n));
    elseif floor(log10(n)) == 2
        data=GetDataSDF(sprintf('0%d.sdf',n));
    end
   n
   data
%% density part
   %% electron
%    n = n+1;
   num_den_ele = data.Derived.Number_Density.el.data;
%    num_den_inj = data.Derived.Number_Density.seventh.data;

   szc = size(num_den_ele);
   %num_den_ele_xy = num_den_ele(:,:,ceil(szc(3)/2));
   %num_den_ele_xz = num_den_ele(:,ceil(szc(2)/2),:);
   %num_den_ele_yz = num_den_ele(ceil(szc(1)/2),:,:);
   %num_den_ele_xy = squeeze(num_den_ele_xy);
   %num_den_ele_xz = squeeze(num_den_ele_xz);
   %num_den_ele_yz = squeeze(num_den_ele_yz);
   x1=data.Derived.Number_Density.el.grid.x;
   y1=data.Derived.Number_Density.el.grid.y;
%    x2=data.Derived.Charge_Density.proton.grid.x;
%    y2=data.Derived.Charge_Density.proton.grid.y;

   %z1=data.Derived.Charge_Density.electron.grid.z;
   figure
%    subplot(2,4,1)
   image(x1(1:end-1),y1(1:end-1),num_den_ele','CDatamapping','scaled');
   % image(num_den','CDatamapping','scaled');
%    image(x1(1:end-1),y1(1:end-1),num_den_inj','CDatamapping','scaled');
   colormap(jet);
   savefig(sprintf('den_ele_xy%d.fig',n))
   
%    image(x2(1:end-1),y2(1:end-1),num_den_pro','CDatamapping','scaled');
   % image(num_den','CDatamapping','scaled');
%    colormap(jet);
%    savefig(sprintf('den_pro_xy%d.fig',n))
   
  % image(x1(1:end-1),z1(1:end-1),num_den_ele_xz','CDatamapping','scaled');
  % colormap(jet);
  % savefig(sprintf('den_ele_xz%d.fig',n))
  % 
  % image(y1(1:end-1),z1(1:end-1),num_den_ele_yz','CDatamapping','scaled');
  % colormap(jet);
  % savefig(sprintf('den_ele_yz%d.fig',n))
%   print('-dbitmap',strcat('den_e',num2str(n),'.bmp'));
%    num_den_beam = data.Derived.Charge_Density.Beam.data;
%    szc = size(num_den_beam);
%    %num_den_ele = num_den_ele(:,ceil(szc(2)/2));
%    num_den_beam = squeeze(num_den_beam);
%    x2=data.Derived.Charge_Density.Beam.grid.x;
%    y2=data.Derived.Charge_Density.Beam.grid.y;
%    subplot(2,4,2)
%    image(x2(1:end-1),y2(1:end-1),num_den_beam','CDatamapping','scaled');
%    % image(num_den','CDatamapping','scaled');
%    colormap(jet);
% %    %set(gca,'Clim',[0 2e24])
%    print('-dbmp',strcat('ele',num2str(n),'.bmp'));
%    %figure
%    %plot(c(:,20))
%    %% line density
%    x_min = 200;
%    x_max = 600;
%    y_min = 76;
%    y_max = 85;
%    num_den_ele_1 = num_den_ele(x_min:x_max,y_min:y_max);
%    figure(n)
%    subplot(2,4,2)
%    image(num_den_ele_1','CDatamapping','scaled');
   e = 1.602176462e-19;
   c = 2.99792458e8;
   epsilon_0 = 8.854187817e-12;
   miu_0     = 1.2566370614e-6;
   me = 9.10938188e-31;
   lamb0 = 0.8e-6;
   E0 = 2*pi*me*c^2/e/lamb0;
%    np = 5e24;             % plasma density
%    wp = e*sqrt(np/epsilon_0/me);
%    l_q0 = c./wp;
%    dense_line = 0;
%    i = 0;
%    for j = y_min:y_max
%        i = i+1;
%        dense_line = num_den_ele_1(:,i)*y1(j)*deta_y+dense_line;   %sum(num_den_ele_1(:,)*y1(100:110));
%    end
%    dense_line = dense_line/np/l_0^2;
%    figure(n)
%    subplot(2,4,3)
%    plot(x1(x_min:x_max),dense_line)
   
   %% proton
  % num_den_pro = data.Derived.Number_Density.proton.data;
  % szc = size(num_den_pro);
  % % b = num_den(:,:,ceil(szc(3)/2));
  % % c = squeeze(b);
  % x2=data.Derived.Number_Density.proton.grid.x;
  % y2=data.Derived.Number_Density.proton.grid.y;
  % figure
  % image(x2(1:end-1),y2(1:end-1),num_den_pro','CDatamapping','scaled');
  % % image(num_den','CDatamapping','scaled');
  % colormap(jet);
  % set(gca,'Clim',[0 2e25])
  % print('-dbmp',strcat('pro',num2str(n),'.bmp'));
   %% carbon
   %num_den_car = data.Derived.Number_Density.carbon.data;
   %szc = size(num_den_car);
   % b = num_den(:,:,ceil(szc(3)/2));
   % c = squeeze(b);
   %x=data.Derived.Number_Density.carbon.grid.x;
   %y=data.Derived.Number_Density.carbon.grid.y;
   %figure
   %image(x(1:end-1),y(1:end-1),num_den_car','CDatamapping','scaled');
   % image(num_den','CDatamapping','scaled');
   %colormap(jet);
   %set(gca,'Clim',[0 5e26])
   %print('-dbmp',strcat('car',num2str(n),'.bmp'));
%% field part
Ex = data.Electric_Field.Ex.data;
Ey = data.Electric_Field.Ey.data;
% Ez = data.Electric_Field.Ez.data;
x = data.Electric_Field.Ex.grid.x;
y = data.Electric_Field.Ex.grid.y;
%z = data.Electric_Field.Ex.grid.z;
sze = size(Ex);
% Ex_xy = squeeze(Ex(:,:,ceil(sze(3)/2)));
% Ex_xz = squeeze(Ex(:,ceil(sze(2)/2),:));
% Ey_xy = squeeze(Ey(:,:,ceil(sze(3)/2)));
% Ey_xz = squeeze(Ey(:,ceil(sze(2)/2),:));
% Ez_xy = squeeze(Ez(:,:,ceil(sze(3)/2)));
% Ez_xz = squeeze(Ez(:,ceil(sze(2)/2),:));
% Ex_yz = squeeze(Ex(ceil(sze(1)/2),:,:));
% Ey_yz = squeeze(Ey(ceil(sze(1)/2),:,:));
% Ez_yz = squeeze(Ez(ceil(sze(1)/2),:,:));



% Bx = data.Magnetic_Field.Bx.data;
% By = data.Magnetic_Field.By.data;
% Bz = data.Magnetic_Field.Bz.data;
% sze = size(Bx);
% 
% u_E  = 0.5*(epsilon_0*(Ex.^2+Ey.^2+Ez.^2));
% u_B  = 0.5*((Bx.^2+By.^2+Bz.^2)/miu_0);
% u    = u_E + u_B;

% Bx_xy = squeeze(Bx(:,:,ceil(sze(3)/2)));
% Bx_xz = squeeze(Bx(:,ceil(sze(2)/2),:));
% By_xy = squeeze(By(:,:,ceil(sze(3)/2)));
% By_xz = squeeze(By(:,ceil(sze(2)/2),:));
% Bz_xy = squeeze(Bz(:,:,ceil(sze(3)/2)));
% Bz_xz = squeeze(Bz(:,ceil(sze(2)/2),:));
% Bx_yz = squeeze(Bx(ceil(sze(1)/2),:,:));
% By_yz = squeeze(By(ceil(sze(1)/2),:,:));
% Bz_yz = squeeze(Bz(ceil(sze(1)/2),:,:));

% Fy = Ey-Bz;
% Fz_xz = Ez_xz+c*By_xz;
% Fz_xy = Ez_xy+c*By_xy;
% x3 = data.Electric_Field.Ey.Core_centre.grid.x;   %data.Derived.Charge_Density.grid.x;
% y3 = data.Electric_Field.Ey.Core_centre.grid.y;
% z = data.Electric_Field.Ez.grid.z;
% figure(n)
% subplot(2,4,3)

% image(x(1:end-1),y(1:end-1),Ex','CDatamapping','scaled');
% % image(Ex','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% savefig(sprintf('u_E%d.fig',n))
% 
% image(x(1:end-1),y(1:end-1),Ey','CDatamapping','scaled');
% % image(Ex','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% savefig(sprintf('u_B%d.fig',n))
% 
% image(x(1:end-1),y(1:end-1),u','CDatamapping','scaled');
% % image(Ex','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% savefig(sprintf('u%d.fig',n))

%% Ex
figure
image(x(1:end-1),y(1:end-1),Ex','CDatamapping','scaled');
% image(Ex','CDatamapping','scaled');
colormap(jet);
colorbar('location','EastOutside')
savefig(sprintf('Ex%d.fig',n))
%print('-dbitmap',strcat('Ex',num2str(n),'.bmp'));

figure
image(x(1:end-1),y(1:end-1),Ey','CDatamapping','scaled');
colormap(jet);
colorbar('location','EastOutside')
%set(gca,'Clim',[-E_m E_m])
savefig(sprintf('Ey%d.fig',n))

% image(x(1:end-1),y(1:end-1),Ez','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% set(gca,'Clim',[-E_m E_m])
% savefig(sprintf('Ez%d.fig',n))

% image(x(1:end-1),y(1:end-1),u','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% %set(gca,'Clim',[-E_m E_m])
% savefig(sprintf('u%d.fig',n))

%% Ey
% subplot(2,4,4)
%image(x(1:end-1),y(1:end-1),Ey_xy','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%% plot(x3(1:end-1),Ex(:,ceil(szc(2)/2)))
%%set(gca,'Clim',[-E_m E_m])
%%print('-dbmp',strcat('ex',num2str(n),'.bmp'));
%savefig(sprintf('Ey_xy%d.fig',n))
%%print('-dbitmap',strcat('Ey',num2str(n),'.bmp'));
%
%image(x(1:end-1),z(1:end-1),Ey_xz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-E_m E_m])
%savefig(sprintf('Ey_xz%d.fig',n))
%
%image(y(1:end-1),z(1:end-1),Ey_yz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-E_m E_m])
%savefig(sprintf('Ey_yz%d.fig',n))
%
%
%%% Ez
%% subplot(2,4,5)
%image(x(1:end-1),y(1:end-1),Ez_xy','CDatamapping','scaled');
%%set(gca,'Clim',[-E_m E_m])
%% image(x3(1:end-1),y3(1:end-1),Ey','CDatamapping','scaled');
%%%% image(Ey','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%savefig(sprintf('Ez_xy%d.fig',n))
%%print('-dbitmap',strcat('Ez',num2str(n),'.bmp'));
%
%image(x(1:end-1),z(1:end-1),Ez_xz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-E_m E_m])
%savefig(sprintf('Ez_xz%d.fig',n))
%
%image(y(1:end-1),z(1:end-1),Ez_yz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-E_m E_m])
%savefig(sprintf('Ez_yz%d.fig',n))

% %% Bx
% % subplot(2,4,5)
% image(x(1:end-1),y(1:end-1),Bx','CDatamapping','scaled');
% % image(x3(1:end-1),y3(1:end-1),Ez','CDatamapping','scaled');
% %%image(Ez','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% savefig(sprintf('Bx%d.fig',n))
% 
% image(x(1:end-1),y(1:end-1),By','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% %set(gca,'Clim',[-2e11 2e11])
% savefig(sprintf('By%d.fig',n))
% 
% image(y(1:end-1),y(1:end-1),Bz','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')
% %set(gca,'Clim',[-2e11 2e11])
% savefig(sprintf('Bz%d.fig',n))

%print('-dbitmap',strcat('Bx',num2str(n),'.bmp'));
% % % set(gca,'Clim',[-1e13 1e13])
% %%% print('-dbmp',strcat('ez',num2str(n),'.bmp'));
% % figure
% % plot(x(1:end-1),c(:,ceil(sze(2)/2)))

%% By 
% subplot(2,4,6)
%image(x(1:end-1),y(1:end-1),By_xy','CDatamapping','scaled');
%% image(Ex','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%savefig(sprintf('By_xy%d.fig',n))
%%print('-dbitmap',strcat('By',num2str(n),'.bmp'));
%
%image(x(1:end-1),z(1:end-1),By_xz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-2e11 2e11])
%savefig(sprintf('By_xz%d.fig',n))
%
%image(y(1:end-1),z(1:end-1),By_yz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-2e11 2e11])
%savefig(sprintf('By_yz%d.fig',n))
%
%%% Bz
%% subplot(2,4,7)
%image(x(1:end-1),y(1:end-1),Bz_xy','CDatamapping','scaled');
%% image(Ex','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%savefig(sprintf('Bz_xy%d.fig',n))
%%print('-dbitmap',strcat('Bz',num2str(n),'.bmp'));
%
%image(x(1:end-1),z(1:end-1),Bz_xz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-2e11 2e11])
%savefig(sprintf('Bz_xz%d.fig',n))
%
%image(y(1:end-1),z(1:end-1),Bx_yz','CDatamapping','scaled');
%colormap(jet);
%colorbar('location','EastOutside')
%%set(gca,'Clim',[-2e11 2e11])
%savefig(sprintf('Bz_yz%d.fig',n))
clear E*
clear B*
clear u*
% subplot(2,4,8)
% image(x(1:end-1),y(1:end-1),Fy','CDatamapping','scaled');
% % image(Ex','CDatamapping','scaled');
% colormap(jet);
% colorbar('location','EastOutside')

% subplot(2,5,8)
% image(x(1:end-1),y(1:end-1),Fz_xz','CDatamapping','scaled');
% % image(Ex','CDatamapping','scaled');
% colormap(jet);

% subplot(2,5,9)
% plot(x(1:end-1),Fz_xy(:,ceil(sze(3)/2)));

%% magnitic field
% Bx = data.Magnetic_Field.Bx.data;
%By = data.Magnetic_Field.By.data;
% Bz = data.Magnetic_Field.Bz.data;
% sze = size(By);
% % b = Ex(:,:,ceil(sze(3)/2));
% % c = squeeze(b);
% x4 = data.Magnetic_Field.By.grid.x;   %data.Derived.Charge_Density.grid.x;
% y4 = data.Magnetic_Field.By.grid.y;
% % z = data.Electric_Field.Ez.grid.z;
% % figure
% % image(x(1:end-1),y(1:end-1),Bx','CDatamapping','scaled');
% % % image(Ex','CDatamapping','scaled');
% % % colormap(jet);
% % % print('-dbmp',strcat('ex',num2str(n),'.bmp'));
% figure
% image(x4(1:end-1),y4(1:end-1),By','CDatamapping','scaled');
% % % image(Ey','CDatamapping','scaled');
% colormap(jet);
% figure
% image(x(1:end-1),y(1:end-1),Bz','CDatamapping','scaled');
% %image(Ez','CDatamapping','scaled');
% colormap(jet);
% % set(gca,'Clim',[-1e13 1e13])
% print('-dbmp',strcat('bz',num2str(n),'.bmp'));

%% particle part
% Px_beam = data.Particles.Px.subset_background.Beam.data;
 try
    a = data.Particles.Px.subset_injection.seventh.data;
 catch
    warning('no particles output');
    a = 0;
 end
sza = size(a)
if sza(1) == 1
   no_par = n
else
   Px_ele  = data.Particles.Px.subset_injection.seventh.data;
   Py_ele  = data.Particles.Py.subset_injection.seventh.data;
%    Pz_ele  = data.Particles.Pz.subset_background.electron.data;
   me = 9.1e-31;
   c = 2.99792458e8;
%c = 3e8;
% Px_beam = Px_beam/me/c;
   Px_ele  = Px_ele/me/c;
   Py_ele  = Py_ele/me/c;
%    Pz_ele  = Pz_ele/mRe/c;

 % gamma = 1./sqrt(1-(vx./c).^2);
 % px = gamma.*gamma/c;   %normalized to mec
 %x5 = data.Particles.Px.subset_background.Beam.grid.x;   %data.Derived.Charge_Density.grid.x;
 %y5 = data.Particles.Px.subset_background.Beam.grid.y;
   x6 = data.Particles.Px.subset_injection.seventh.grid.x;   %data.Derived.Charge_Density.grid.x;
   y6 = data.Particles.Px.subset_injection.seventh.grid.y;
%   z6 = data.Particles.Px.subset_background.electron.grid.z;
   gamma=data.Particles.Gamma.subset_injection.seventh.data;

 % z = data.Particles.Gamma.subset_background.electron.grid.z;
 %i  = i+1;
 %subplot(2,4,6)
 %plot(x5,Px_beam,'.','markersize',0.1);
   subplot(3,3,1)
   plot(x6,Px_ele,'.','markersize',0.1);
   subplot(3,3,2)
   plot(y6,Px_ele,'.','markersize',0.1);
   subplot(3,3,3)
   plot(x6,y6,'.','markersize',0.1);
   subplot(3,3,4)
   plot(x6,gamma,'.','markersize',0.1);
%    subplot(3,3,5)
figure
   plot(y6,Py_ele,'.','markersize',0.1);
   %subplot(3,3,6)
   %plot(z6,Py_ele,'.','markersize',0.1);
%    subplot(3,3,7)
%    plot(x6,Pz_ele,'.','markersize',0.1);
%    subplot(3,3,8)
%    plot(y6,Pz_ele,'.','markersize',0.1);
   %subplot(3,3,9)
   %plot(z6,Pz_ele,'.','markersize',0.1);
   savefig(sprintf('px-py-pz%d.fig',n))
end
 %ylim([-12.0e-20,12e-20]);
 %print('-dbmp',strcat('px',num2str(n),'.bmp'));
% clear data
%clear a
 %% particle part
end
% matlabpool close
toc
