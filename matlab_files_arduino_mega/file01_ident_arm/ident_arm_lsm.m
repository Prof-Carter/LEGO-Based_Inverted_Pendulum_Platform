
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
Tf1 = 0.04;

n = length(t);

% ==================
dtheta1(1) = 0;
for i = 2:n-1
    dtheta1(i) = (theta1(i+1) - theta1(i-1))/(2*h);
end
dtheta1(n) = dtheta1(n-1);

dtheta1 = dtheta1';

% ==================
ddtheta1(1) = 0;
for i = 2:n-1
    ddtheta1(i) = (dtheta1(i+1) - dtheta1(i-1))/(2*h);
end
ddtheta1(n) = ddtheta1(n-1);

ddtheta1 = ddtheta1';

% ----------------------------------------------
t_begin = 1;
t_end   = 5;

data_begin = t_begin/h + 1;
data_end   = t_end/h   + 1;

t_data = t(data_begin:data_end) - t_begin;

  theta1_data =   theta1(data_begin:data_end);
 dtheta1_data =  dtheta1(data_begin:data_end);
ddtheta1_data = ddtheta1(data_begin:data_end);

r_data = r(data_begin:data_end);
u_data = u(data_begin:data_end);


% ----------------------------------------------
clear t theta1 r u

t = t_data;

  theta1 =   theta1_data;
 dtheta1 =  dtheta1_data;
ddtheta1 = ddtheta1_data;

r = r_data;
u = u_data;

rc = r(end);

% ----------------------------------------------
M1 = ddtheta1;
M2 = dtheta1;
M3 = sign(dtheta1);
M  = [ M1  M2  M3 ];
N  = u;

sysGf1 = tf(1,[Tf1 1])^3;

Mf1 = lsim(sysGf1,M1,t);
Mf2 = lsim(sysGf1,M2,t);
Mf3 = lsim(sysGf1,M3,t);
Mf  = [Mf1 Mf2 Mf3];
Nf  = lsim(sysGf1,N,t);

p = inv(Mf'*Mf)*Mf'*Nf;

alpha1 = p(1);
mu1    = p(2);
delta1 = p(3);

fprintf('alpha1 = %3.2e;\n',alpha1)
fprintf('mu1    = %3.2e;\n',mu1)
fprintf('delta1 = %3.2e;\n',delta1)

% ----------------------------------------------
sim('sim_arm_pcont')

figure(1)
set(gcf,'Color','white','PaperType','A3')
set(gcf,'Position',[100 100 900 360]) % [x0 y0 width height]
subplot('Position',[0.15 0.2 0.8 0.7]) % [left bottom width height]

plot([0 3],[rc rc]*180/pi,'k')
hold on
p2 = plot(t,theta1*180/pi,'Color',myred,'LineWidth',1.5);
p1 = plot(t_sim,theta1_sim*180/pi,'--','Color',mygreen,'LineWidth',1.5);
hold off

set(gca,'FontName','arial','FontSize',20)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',22)
ylabel('$${\theta}_{1}$$ [deg]', 'interpreter', 'latex','FontSize',22)

xlim([0 2])
ylim([0 160])

set(gca,'XTick',0:0.2:2)
set(gca,'YTick',0:30:180)

legend([p2 p1],{'Experiment','Nonlinear simulation'},'Location','southeast')
set(legend, 'interpreter','latex','FontSize',22)

grid on
