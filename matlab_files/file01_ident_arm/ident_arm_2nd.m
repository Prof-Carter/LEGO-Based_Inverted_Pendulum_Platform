clear
format compact

myred   = [232   71   70]/255;
myblue	= [  0  112  192]/255;
mygreen	= [ 12  162  116]/255;
mypink	= [255  102  153]/255;
myskyblue    = [  0  176  240]/255;
mydarkpink   = [227   45  145]/255;
mydarkyellow = [228  131   18]/255;

% ----------------------------------------------
load arm_data

n = length(t);

% ----------------------------------------------
t_begin = 1;
t_end   = 5;

data_begin = t_begin/h + 1;
data_end   = t_end/h   + 1;

t_data      = t(data_begin:data_end) - t_begin;
theta1_data = theta1(data_begin:data_end);
r_data      = r(data_begin:data_end);
u_data      = u(data_begin:data_end);

% ----------------------------------------------
clear t u theta1 r

t      = t_data;
theta1 = theta1_data;
r      = r_data;

rc = r(end);

% ----------------------------------------------
% ++++++++++ Step 1 ++++++++++
[theta1_max i_max] = max(theta1);
Amax = theta1_max - rc;
Tp   = t(i_max);

fprintf('Tp   = %3.2e;\n',Tp)
fprintf('Amax = %3.2e;\n',Amax)

% ++++++++++ Step 2 ++++++++++
xi   = - (1/Tp)*log(Amax/rc);
wn   = sqrt((pi/Tp)^2 + xi^2);      % 固有角周波数
zeta = xi/wn;                       % 減衰係数

fprintf('wn   = %3.2e;\n',wn)
fprintf('zeta = %3.2e;\n',zeta)

% ++++++++++ Step 3 ++++++++++
alpha1 = kP/wn^2;                   % alpha1 の同定
mu1    = 2*zeta*wn*alpha1;          % mu1 の同定

fprintf('alpha1 = %3.2e;\n',alpha1)
fprintf('mu1    = %3.2e;\n',mu1)
fprintf('delta1 = 0;\n')

% -----------------------
t_sim = 0:0.001:3;                  % 同定されたパラメータを用いてシミュレーション
theta1_sim = rc*step(wn^2,[1 2*zeta*wn wn^2],t_sim);

figure(1)
movegui('north')
set(gcf,'Color','white','PaperType','A3')
set(gcf,'Position',[100 400 900 360]) % [x0 y0 width height]
subplot('Position',[0.15 0.2 0.8 0.7]) % [left bottom width height]

plot([Tp Tp],[0 theta1_max]*180/pi,'k')
hold on
plot([0 Tp],[theta1_max theta1_max]*180/pi,'k')
plot([0 3],[rc rc]*180/pi,'k')

p2 = plot(t,theta1*180/pi,'Color',myred,'LineWidth',1.5);
p1 = plot(t_sim,theta1_sim*180/pi,'--','Color',mygreen,'LineWidth',1.5);

plot(Tp,theta1_max*180/pi,'o','LineWidth',1.5,'MarkerSize',8,'Color',myblue)
hold off

set(gca,'FontName','arial','FontSize',20)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',22)
ylabel('$${\theta}_{1}$$ [deg]', 'interpreter', 'latex','FontSize',22)

xlim([0 2])
ylim([0 160])

set(gca,'XTick',0:0.2:2)
set(gca,'YTick',0:30:180)

legend([p2 p1],{'Experiment','Linear simulation'},'Location','southeast')
set(legend, 'interpreter','latex','FontSize',22)

grid on