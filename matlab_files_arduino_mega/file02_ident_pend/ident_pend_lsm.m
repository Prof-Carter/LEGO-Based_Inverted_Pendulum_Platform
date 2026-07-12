clear
format compact
close all

myred   = [232   71   70]/255;
myblue	= [  0  112  192]/255;
mygreen	= [ 12  162  116]/255;
mypink	= [255  102  153]/255;
myskyblue    = [  0  176  240]/255;
mydarkpink   = [227   45  145]/255;
mydarkyellow = [228  131   18]/255;

% -----------------------
g  = 9.81e+00;
Tf2 = 0.05;
% -----------------------
load pend_data

tmin = 4;
tmax = 14;

% -----------------------
figure(1)
set(gcf,'Color','white','PaperType','A3')
set(gcf,'Position',[100 100 900 540]) % [x0 y0 width height]
% set(gcf,'Position',[100 400 900 360]) % [x0 y0 width height]
subplot('Position',[0.15 0.2 0.8 0.7]) % [left bottom width height]

area([0 tmin], 110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0)
hold on
area([0 tmin], -110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0)
area([tmax 15],110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0)
area([tmax 15],-110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0)
plot([tmin tmin],[-110 110],'Color',mydarkyellow)

q1 = stairs(t,phi2*180/pi,'Color',myred,'LineWidth',1.5);
hold off

xlim([0 tmax+1])
ylim([-110 110])

set(gca,'XTick',0:1:tmax+1)
set(gca,'YTick',-90:30:90)

set(gca,'FontName','arial','FontSize',20)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',22)
ylabel('$${\phi}_{2}$$ [deg]', 'interpreter', 'latex','FontSize',22)

grid on

% =========================================
k = length(t);
% -----------------------
dphi2(1) = (- 3*phi2(1) + 4*phi2(2) - phi2(3))/(2*h);
for i = 2:k-1
    dphi2(i) = (phi2(i+1) - phi2(i-1))/(2*h);
end
dphi2(k) = (phi2(k-2) - 4*phi2(k-1) + 3*phi2(k))/(2*h);
% -----------------------
ddphi2(1) = (- 3*dphi2(1) + 4*dphi2(2) - dphi2(3))/(2*h);
for i = 2:k-1
    ddphi2(i) = (dphi2(i+1) - dphi2(i-1))/(2*h);
end
ddphi2(k) = (dphi2(k-2) - 4*dphi2(k-1) + 3*dphi2(k))/(2*h);
% -----------------------

% =========================================
M1 = ddphi2;
M2 = dphi2;
N  = - g*sin(phi2);

sysGf2 = tf(1,[Tf2 1])^3;

Mf1 = lsim(sysGf2,M1,t);
Mf2 = lsim(sysGf2,M2,t);
Nf  = lsim(sysGf2,N ,t);

% =========================================
t_data = t;

phi2_data   = phi2;
dphi2_data  = dphi2;
ddphi2_data = ddphi2;

M1_data = M1;
M2_data = M2;
N_data  = N;

Mf1_data = Mf1;
Mf2_data = Mf2;
Nf_data  = Nf;

clear t phi2 dphi2 ddphi2
clear M1 M2 N Mf1 Mf2 Nf
% -----------------------
k = 0;
for i = 1:length(t_data)
    if t_data(i) >= tmin & t_data(i) <= tmax
        k = k + 1;
        t(k,1)   = t_data(i) - tmin;

        phi2(k,1)   = phi2_data(i);
        dphi2(k,1)  = dphi2_data(i);
        ddphi2(k,1) = ddphi2_data(i);

        Mf1(k,1) = Mf1_data(i);
        Mf2(k,1) = Mf2_data(i);
        Nf(k,1)  = Nf_data(i);

        M1(k,1) = M1_data(i);
        M2(k,1) = M2_data(i);
        N(k,1)  = N_data(i);
    end
end

% =========================================
Mf = [ Mf1  Mf2 ];

p2 = inv(Mf'*Mf)*Mf'*Nf;

alpha2 = p2(1);
mu2    = p2(2);

fprintf('alpha2 = %3.2e;\n',alpha2)
fprintf('mu2    = %3.2e;\n',mu2 )

% -----------------------
[phi20 i] = max(phi2);
t0 = t(i);

sim('sim_pend_free')

% ================================
figure(1)
hold on
q2 = plot(tsim+tmin,phi2sim*180/pi,'--','Color',mygreen,'LineWidth',1.5);
plot(tsim(1)+tmin,phi2sim(1)*180/pi,'o','Color',myblue,'LineWidth',1.5,'MarkerSize',8)
hold off

legend([q1 q2],{'Experiment$\ \ $','Nonlinear simulation'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',22)
