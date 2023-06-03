clc
clear
close all
e = 1.602176462e-19;
c = 2.99792458e8;
epsilon_0 = 8.854187817e-12;
miu_0     = 1.2566370614e-6;
me = 9.10938188e-31;

%% 判断挑选粒子的范围
n=150;
if floor(log10(n)) == 0
        refer=GetDataSDF(sprintf('000%d.sdf',n));
    elseif floor(log10(n)) == 1
        refer=GetDataSDF(sprintf('00%d.sdf',n));
    elseif floor(log10(n)) == 2
        refer=GetDataSDF(sprintf('0%d.sdf',n));
 end
   n
num_den_inj = refer.Derived.Number_Density.ionisation.data;
z2 = flipud(rot90(num_den_inj))/1e6;
x_1=refer.Derived.Number_Density.ionisation.grid.x;
y_1=refer.Derived.Number_Density.ionisation.grid.y;
figure
image(x_1(1:end-1,:),y_1(1:end-1,:),z2,'CDatamapping','scaled');%通过在创建图像时将 CDataMapping 属性设置为 'scaled'，将值的范围缩放到当前颜色图的完整范围。
axis xy
 load('myclmp.mat')
 colormap(CustomColormap1)
 colorbar 
 xlabel('x(m)','FontSize',15)
 ylabel('y(m)','FontSize',15)
 set(gca,'XLim',[2.98 2.985]*1e-3)
 set(gca,'YLim',[-1 1]*1e-5)
 x_2=refer.Particles.ID.subset_injection.ionisation.grid.x;
 index=find(x_2>2.981e-3 & x_2<2.984e-3);
 hold on
 plot([2.981e-3 2.981e-3],[-1e-5,0],'r-','LineWidth',2)
 plot([2.984e-3 2.984e-3],[-1e-5,0],'r-','LineWidth',2)
 hold off

 %% 随机不重复的从数组中选取指定个数ID
 Id_r1 = refer.Particles.ID.subset_injection.ionisation.data(index);
 Id_r1=Id_r1(randperm(numel(Id_r1),2000));%这里的2000表示选取的粒子个数
 for n =15
    if floor(log10(n)) == 0
        data=GetDataSDF(sprintf('000%d.sdf',n));
    elseif floor(log10(n)) == 1
        data=GetDataSDF(sprintf('00%d.sdf',n));
    elseif floor(log10(n)) == 2
        data=GetDataSDF(sprintf('0%d.sdf',n));
    end
    n
   Idcycle = data.Particles.ID.subset_injection.ionisation.data;
   IDkkk=intersect(Id_r1,Idcycle);  
   Id_r1=IDkkk;
 end

 %% 导出粒子信息(x,y,z,βx,βy,βz,t,γ)，存为track_*.dat的文件
 for n = 15:1:150
    if floor(log10(n)) == 0
        data=GetDataSDF(sprintf('000%d.sdf',n));
    elseif floor(log10(n)) == 1
        data=GetDataSDF(sprintf('00%d.sdf',n));
    elseif floor(log10(n)) == 2
        data=GetDataSDF(sprintf('0%d.sdf',n));
    end
   data;
   n
   Id1 = data.Particles.ID.subset_injection.ionisation.data(:,1);
   [Idreal, ir, ib] = intersect(Id_r1,Id1);
   x1 = data.Particles.ID.subset_injection.ionisation.grid.x(:,1);
   y1 = data.Particles.ID.subset_injection.ionisation.grid.y(:,1);
   z1=data.Particles.ID.subset_injection.ionisation.grid.y(:,1);
   z1=z1.*0;
   Px_ele  = data.Particles.Px.subset_injection.ionisation.data(:,1);
   Py_ele  = data.Particles.Py.subset_injection.ionisation.data(:,1);
   Pz_ele  = data.Particles.Pz.subset_injection.ionisation.data(:,1);
   Gamma=data.Particles.Gamma.subset_injection.ionisation.data(:,1);
   T = data.time;
   format long
   Px  = Px_ele./Gamma./me./c;
   Py  = Py_ele./Gamma./me./c;
   Pz  = Pz_ele./Gamma./me./c;
   for j = 1:size(Idreal)            
   Ytrack(n,j) = y1(ib(j));
   Xtrack(n,j) = x1(ib(j));
   Ztrack(n,j) = z1(ib(j));
   Pxtrack(n,j) = Px(ib(j));
   Pytrack(n,j) = Py(ib(j));
   Pztrack(n,j) = Pz(ib(j));
   Gammatrack(n,j) = Gamma(ib(j));
   time(n,j) = T;
   end
 end
 format long e
 for i= 1:size(Idreal) 
 trace=[Ytrack(15:1:150,i) Ztrack(15:1:150,i) Xtrack(15:1:150,i) Pytrack(15:1:150,i) Pztrack(15:1:150,i) Pxtrack(15:1:150,i) time(15:1:150,i) Gammatrack(15:1:150,i)];
 str= ['track_' num2str(i) '.dat'];
 fid=fopen(str,'w');
 [m,n2] = size(trace);
   for l= 1:1:m
       for t = 1:1:n2
           if t == n2
               fprintf(fid,'%0.15e\n',trace(l,t));
           else
               fprintf(fid,'%0.15e\t',trace(l,t));
           end
       end
   end
 fclose(fid);
 end
%% 粒子的轨迹
figure
c=colormap(lines(i));
j=1;
sp=1;
for q=1:i
    str=['track_' num2str(q) '.dat'];
    x_1=importdata(str);
    x=x_1(:,3);
    y=x_1(:,1);
    i_n=find(x>1*1e-3);
    y_max=max(y(i_n));
    if(y_max>1*1e-6)%挑选指定振幅的粒子的轨迹
    plot(x*1e3,y*1e6,'color',c(q,:))
    p(1,sp)=q;
    sp=sp+1;
    end
    j=j+1;
    hold on
end
box on
ylabel('y(μm)')
xlabel('x(mm)')
set(gca,'XLim',[0.5,3],'XTick',[0.5:0.5:3])
set(gca,'YLim',[-4,4],'YTick',[-4:1:4])
set(gca,'fontsize',15,'fontname','Times New Roman')
%% 保存被挑选出的特定振幅粒子的dat数据
for ii=1:length(p)
    str2=['track_' num2str(p(ii)) '.dat'];
    x_2=importdata(str2);
    trace2=x_2;
    str= ['spe_track' num2str(ii) '.dat'];
    fid=fopen(str,'w');
 [m,n2] = size(trace2);
   for l= 1:1:m
       for t = 1:1:n2
           if t == n2
               fprintf(fid,'%0.15e\n',trace2(l,t));
           else
               fprintf(fid,'%0.15e\t',trace2(l,t));
           end
       end
   end
   fclose(fid);
end