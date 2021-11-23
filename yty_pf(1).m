%开门三件事：清变量、关窗口、清命令行
clear;
close all;
clc;

%变量初始化
n=100;               %粒子数
step=100;            %步数
delta_t=0.00001;     %时间微分
x=zeros(100,1);      %x坐标真实值
y=zeros(100,1);      %y坐标真实值
theta=zeros(100,1);  %方向角
yawrate=zeros(100,1);%方向角加速度
for i=1:100
    theta(i)=(0.01*(50-i))^3;
end
for  i=2:100
    yawrate(i)=(theta(i)-theta(i-1))/delta_t;
end
x(1)=0.1;            %x坐标初始值
y(1)=0.2;            %y坐标初始值
velocity=0.1:0.1:10;         %速度大小
Q=1;                         %过程噪声方差
R=1;                         %观测噪声方差
v=sqrt(R)*randn(step,1);     %观测噪声

%生成真实数据
for i=2:100
    if abs(yawrate(i))<0.001
        theta(i)=theta(i-1);
        x(i)=x(i-1)+velocity(i)*cos(theta(i-1))*delta_t;
        y(i)=y(i-1)+velocity(i)*sin(theta(i-1))*delta_t;
    else
        theta(i)=theta(i-1)+yawrate(i)*delta_t;
        x(i)=x(i-1)+(velocity(i)/yawrate(i))*(sin((theta(i-1)+(yawrate(i)*delta_t))-sin(theta(i-1))));
        y(i)=y(i-1)+(velocity(i)/yawrate(i))*(-cos((theta(i-1)+(yawrate(i)*delta_t))+cos(theta(i-1))));
    end
end

%生成观测数据
zx=x+v; %观测值=真实值+当前步观测噪声
zy=y+v; %观测值=真实值+当前步观测噪声

%设置x粒子集合
xpar=zeros(n,step);  %粒子矩阵
xpf=zeros(n,step);   %滤波后的粒子
xw=zeros(n,step);    %粒子权重
%注：由于观测方程为z=x+v，粒子值即观测值，故不再设置粒子观测值数组
xpf(:,1)=x(1)+sqrt(Q)*randn(n,1);

%PF
for k=2:step
    %1.预测（投粒子）
    for i=1:n
        if abs(yawrate(k))<0.001
            xpar(i,k)=xpf(i,k-1)+velocity(k)*cos(theta(k-1))*delta_t;
        else
            xpar(i,k)=xpf(i,k-1)+...
            (velocity(k)/yawrate(k))*(sin((theta(k-1)+(yawrate(k)*delta_t))-sin(theta(k-1))));
        end
    end
    %2.计算权重
    for i=1:n
        xw(i,k)=exp(-.5*R^(-1)*(zx(k,1)-xpar(i,k))^2);
    end
    xw(:,k)=xw(:,k)./sum(xw(:,k)); %权值归一化
    %3.重采样
    outIndex=randomR(xw(:,k));
    xpf(:,k)=xpar(outIndex,k);
end
%4.粒子滤波结果
xmean=mean(xpf);

%设置y粒子集合
ypar=zeros(n,step);  %粒子矩阵
ypf=zeros(n,step);   %滤波后的粒子
yw=zeros(n,step);    %粒子权重
%注：由于观测方程为z=x+v，粒子值即观测值，故不再设置粒子观测值数组
ypf(:,1)=y(1)+sqrt(Q)*randn(n,1);

%PF
for k=2:step
    %1.预测
    for i=1:n
        if abs(yawrate(k))<0.001
            ypar(i,k)=ypf(i,k-1)+velocity(k)*sin(theta(k-1))*delta_t;
        else
            ypar(i,k)=ypf(i,k-1)+...
            (velocity(k)/yawrate(k))*(-cos((theta(k-1)+(yawrate(k)*delta_t))+cos(theta(k-1))));
        end
    end
    %2.计算权重
    for i=1:n
        yw(i,k)=exp(-.5*R^(-1)*(zy(k,1)-ypar(i,k))^2);
    end
    yw(:,k)=yw(:,k)./sum(yw(:,k)); %权值归一化
    %3.重采样
    outIndex=randomR(yw(:,k));
    ypf(:,k)=ypar(outIndex,k);
end
%4.粒子滤波结果
ymean=mean(ypf);

%粒子滤波结果图
hold on;
title('粒子滤波结果图');
xlabel('x坐标');
ylabel('y坐标');
plot(x,y);
% plot(zx,zy);
% plot(xmean,ymean);
legend('真实位置坐标','观测位置坐标','滤波位置坐标');
hold off;

%滤波结果分析
e_pf=zeros(100,1); %滤波值和真实值的位置偏差
e_z =zeros(100,1); %观测值和真实值的位置偏差
sumdis_pf=0;       %滤波值和真实值的位置偏差总和
sumdis_z =0;       %观测值和真实值的位置偏差总和
for i=1:step
    e_pf(i)=sqrt((xmean(i)-x(i))^2+(ymean(i)-y(i))^2);
    e_z(i) =sqrt((zx(i)-x(i))^2+(zy(i)-y(i))^2);
    sumdis_pf=sumdis_pf+e_pf(i);
    sumdis_z =sumdis_z+e_z(i);
end
sigma_pf=var(e_pf);
sigma_z =var(e_z);
fprintf('滤波值和真实值的位置偏差方差：%f\n',sigma_pf);
fprintf('观测值和真实值的位置偏差方差：%f\n',sigma_z);
fprintf('滤波值和真实值的位置偏差总和：%f\n',sumdis_pf);
fprintf('观测值和真实值的位置偏差总和：%f\n',sumdis_z);