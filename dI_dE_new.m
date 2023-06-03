clc
close all
clear
ne=5e18;%立方cm
 me=9.10938215e-31;%kg
 e = 1.602176462e-19;
C = 2.99792458e+8;
h_bar=6.582119514e-16;%单位ev*s
h=4.1356676969e-15;%ev.s
epsilon_0 = 8.854187817e-12;

da=1:1000;%表示将1000个粒子的辐射能求和
for s=1:length(da)
    str=['track_' num2str(da(s)) '.dat'];
    x_1=importdata(str);
    x_g1=x_1(1:end,3);
    y_g1=x_1(1:end,1);
    gamma=x_1(:,8);
    inde=find(x_g1>1e-3);%计算1mm后的辐射情况

    x_g2=x_g1(inde);
    y_g2=y_g1(inde);
    gamma=gamma(inde);
    d=findpeaks(y_g1);
IndMin=find(diff(sign(diff(y_g2)))>0)+1;   %获得局部最小值的位置
IndMax=find(diff(sign(diff(y_g2)))<0)+1;   %获得局部最大值的位置
% plot(x_g1,y_g1,'color','b','linewidth',2);
% hold on;
% box on;
% plot(x_g2(IndMin),y_g2(IndMin),'r+');
% plot(x_g2(IndMax),y_g2(IndMax),'r+');
% xlabel('x(m)','FontSize',15);
% ylabel('y(m)','FontSize',15);
% ylim([-6 6]*1e-6)
for p1=1:length(IndMin)
    q1=IndMin(p1);
    rB1=abs(y_g2(q1));
    rB1=rB1.*1e6;
    gama=gamma(q1);
    Ec=5.24e-24.*gama^2.*ne.*rB1;%kev
    wc=Ec*1e3/h_bar;
    aa=1;
    for E=0:1:100%kev
        w=E*1e3/h_bar;
        y=@(x)(besselk(5/3,x));
        yy=integral(y,w/wc,inf);
        dI_w=2*sqrt(3)*(e^2/C)*gama*w/wc*yy;
    d_E=5*(3*e^2)/(2*pi^3*h_bar*C*epsilon_0)*gama^2*(E/Ec)^2*(besselk(2/3,E/Ec).^2);%这的5表示电子的振荡周期
    n2(p1,aa)=dI_w;
    n3(p1,aa)=d_E;
    aa=aa+1;
    end
end

for p2=1:length(IndMax)
    q2=IndMax(p2);
    rB1=abs(y_g2(q2));
    rB1=rB1.*1e6;
    gama=gamma(q2);
    Ec=5.24e-24.*gama^2.*ne.*rB1;%kev
    wc=Ec*1e3/h_bar;
    aa=1;
    for E=0:1:100%kev
        w=E*1e3/h_bar;
        y=@(x)(besselk(5/3,x));
        yy=integral(y,w/wc,inf);
        dI_w=2*sqrt(3)*(e^2/C)*gama*w/wc*yy;
    d_E=5*(3*e^2)/(2*pi^3*h_bar*C*epsilon_0)*gama^2*(E/Ec)^2*(besselk(2/3,E/Ec).^2);%这的5表示电子的振荡周期
    n2(p2+length(IndMin),aa)=dI_w;
    n3(p2+length(IndMin),aa)=d_E;
    aa=aa+1;
    end
end
num(s,:)=sum(n3);
num2(s,:)=sum(n2);
n3=[];
n2=[];
end
E=0:1:100;
num_1=sum(num);
num_2=sum(num2);
plot(E,num_1,'k-','LineWidth',2)
hold on
aa=find(num_1==max(num_1));
E_max=E(aa);
plot([E_max,E_max],[0,num_1(aa)],'r--')
text(E_max,0,num2str(E_max),'Color','r','FontSize',15)
xlabel('E(kev)','FontSize',15);
ylabel('d^{2}I/dE\Omega','FontSize',15);
set(gca,'fontname','times new roman')
figure
plot(E,num_2,'k-','LineWidth',2)
hold on
aa=find(num_2==max(num_2));
E_max=E(aa);
plot([E_max,E_max],[0,num_2(aa)],'r--')
text(E_max,0,num2str(E_max),'Color','r','FontSize',15)
xlabel('E(kev)','FontSize',15);
ylabel('dI/d\omega','FontSize',15);
set(gca,'fontname','times new roman')