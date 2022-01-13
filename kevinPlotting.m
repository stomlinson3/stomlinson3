clear all
% close all
%ap

%
% x=(-pi/2+0.01):0.001:(pi/2-0.01);
% y=tan(x);
% figure(1)
% clf
% plot(x,y)
% ylabel('tan(x)')
% ylim([-100,100]);
%
% figure(iiiiiiiiiiiiiiiiiiiiiiiiiiiq2)

% clf
% plot(y,atan(y))
% wasd


% % % sw=load('../space_weather_2021-11-02.csv')
% % % tau=8400;
% % % ss=10000;
% % % dt=0.1
% % %
% % % b(1)=0;
% % % t(1)=0;
% % % ttau(1)=tau;
% % %
% % % for i =2:1:10000
% % %     r=randn(1);
% % %
% % %     b0=b(i-1)+((b(i-1) + sqrt(2)*r*ss))*dt/tau;
% % %
% % %
% % %     t(i)=t(i-1)+dt;
% % %
% % %     ttau(i)=((b0-b(i-1))/((b(i-1) + sqrt(2)*r*ss)*dt))^(-1);
% % %     sss(i)=(b0-b(i-1))*tau/dt-b(i-1)/(sqrt(2)*r);
% % %     b(i)=b0;
% % %
% % % end
% % %
% % % figure(1)
% % % clf
% % % plot(t,b)
% % %
% % % figure(2)
% % % clf
% % % plot(t,ttau)
% % %
% % % figure(3)
% % % clf
% % % plot(t,sss)
% % % sw(:,2)=sw(:,2)-sum(sw(:,2))/length(sw(:,2));
% % % sw(:,3)=sw(:,3)-sum(sw(:,3))/length(sw(:,3));
% % % sw(:,1)=sw(:,1)-sw(1,1)
% % % sw0(:,2:3)=interp1(sw(:,1), sw(:,2:3), t');
% % % sw0(:,1)=t';
% % % sw=sw0;
% % % errsw(:,1)=est(16,:)'-sw(:,2);
% % % errsw(:,2)=est(16,:)'-sw(:,3);
% % % sw2=[];
% % % n=1;
% % % for i =1:1:length(sw(:,1))
% % %     if sw(i,1)>min(t) && sw(i,1)<=max(t)
% % %         sw2(n,:)=sw(i,:);
% % %         n=n+1
% % %     end
% % % end
% % %
% % % sw=sw2;
% % % figure(2)
% % % clf
% % % plot(sw(:,1),sw(:,2))
% % % hold on
% % % plot(sw(:,1),sw(:,3))
% % % hold on
% % % plot(t,est(16,:))
% % % legend('sw1','sw2','FOGM Est')
% % % title('FOGM estimate with GPS on: mean removed')
% % %
% % % figure(2)
% % % clf
% % % plot(t',errsw(:,1))
% % % hold on
% % % plot(t',errsw(:,2))
% % % hold on
% % % plot(t,3.*cov(16,:),'r')
% % % hold on
% % % plot(t,-3.*cov(16,:),'r')
% % % title('error btw sw and FOGM estimate: GPS on')
% % % legend('mag1','mag2','3 sigma')

% % % adfdsafds

cov0=load('Cov.dat');
est0=load('Est.dat');
pva0=load('Pva.dat');
err0=load('err.dat');
t0=load('Time.dat');
n0=length(t0);
n1=length(t0)/n0

cov=zeros(16,n0);
err=zeros(16,n0);
est=zeros(16,n0);
pva=zeros(16,n0);
t=zeros(1,n0);

errA=zeros(16,n0);

n2=1;
n3=1;
for i=1:1:n1
    for j =1:1:n0
        for k=1:1:16
            cov(k,j)=cov(k,j)+cov0(n2)/n1;
            
            error=err0(n2);
            
            err(k,j)=err(k,j)+error/n1;
            
            errA(k,j)=errA(k,j)+sqrt(error*error)/n1;
            
            pva(k,j)=pva(k,j)+pva0(n2)/n1;
            est(k,j)=est(k,j)+est0(n2)/n1;
            n2 = n2 + 1;
        end
        t(1,j)=t(1,j) + t0(n3)/n1;
        n3 = n3 + 1;
    end
end

est(16,1)=est(16,2);
% est(16,:)=est(16,:)-sum(est(16,:))/length(est(16,:));
t(1,:)=t(1,:)-t(1,1);

color = 'k'
linewidth = 10.0

figure(1)
clf
n=1;
for j=1:1:3
    for i=1:1:5
        subplot(3,5,n)
        plot(t,err(3*(i-1)+j,:),color)
        hold on
        plot(t,cov(3*(i-1)+j,:),color)
        hold on
        plot(t,-cov(3*(i-1)+j,:),color)
        hold on
        set(gca,'Fontsize',16)
        n=n+1;
    end
end

subplot(3,5,1)
title('position')
subplot(3,5,2)
title('velocity')
subplot(3,5,3)
title('attitude')
subplot(3,5,4)
title('accel bias')
subplot(3,5,5)
title('gyro bias')
subplot(3,5,1)
ylabel('x')
subplot(3,5,6)
ylabel('y')
subplot(3,5,11)
ylabel('z')
xlabel('time(s)')

figure(2)
clf
plot(pva(2,:)*180.0/pi,pva(1,:)*180.0/pi,'k')
hold on
plot(est(2,:)*180.0/pi,est(1,:)*180.0/pi,'r')


figure(3)
clf
n=1;
for j=1:1:3
    for i=1:1:5
        subplot(3,5,n)
        plot(t,errA(3*(i-1)+j,:),color)
        hold on
        plot(t,cov(3*(i-1)+j,:)./3.0,color)
        hold on
        set(gca,'Fontsize',16)
        n=n+1;
    end
end

subplot(3,5,1)
title('position')
subplot(3,5,2)
title('velocity')
subplot(3,5,3)
title('attitude')
subplot(3,5,4)
title('accel bias')
subplot(3,5,5)
title('gyro bias')
subplot(3,5,1)
ylabel('x')
subplot(3,5,6)
ylabel('y')
subplot(3,5,11)
ylabel('z')
xlabel('time(s)')

figure(4)
clf
subplot(2,1,1)

plot(t,err(1,:),color)
hold on
plot(t,cov(1,:),color)
hold on
plot(t,-cov(1,:),color)
hold on
title('errors (m)')
ylabel('north')
set(gca,'Fontsize',16)

subplot(2,1,2)
plot(t,err(2,:),color)
hold on
plot(t,cov(2,:),color)
hold on
plot(t,-cov(2,:),color)
hold on
ylabel('east')
xlabel('time (s)')
set(gca,'Fontsize',16)









est(1:2,:)=est(1:2,:).*180/pi;
pva(1:2,:)=pva(1:2,:).*180/pi;
figure(5)
clf
n=1;
for j=1:1:3
    for i=1:1:3
        subplot(3,3,n)
        plot(t,est(3*(i-1)+j,:),'r')
        hold on
        plot(t,pva(3*(i-1)+j,:),'b')
        hold on
        set(gca,'Fontsize',16)
        n=n+1;
    end
end


subplot(3,3,1)
title('position')
legend('updated','truth')
subplot(3,3,2)
title('velocity')
subplot(3,3,3)
title('attitude')
subplot(3,3,1)
ylabel('x')
subplot(3,3,4)
ylabel('y')
subplot(3,3,7)
ylabel('z')
xlabel('time(s)')
set(gca,'Fontsize',16)






set(gca,'Fontsize',16)

figure(6)
clf
plot(t,err(16,:),'b')
hold on
plot(t,cov(16,:),'b')
hold on
plot(t,-cov(16,:),'b')
hold on
title('error mag bias')
ylabel('nT')
xlabel('time(s)')
set(gca,'Fontsize',16)



figure(7)
clf
plot(t(2:length(est(1,:))),est(16,2:length(est(1,:))),'r')
hold on
plot(t,pva(16,:),'b')
hold on
title('mag bias')
legend('estimate','true')
set(gca,'Fontsize',16)

avg=[];
for i =1:1:16
    avg(i)=sum(err(i,:))/length(err(i,:));
end
avg

