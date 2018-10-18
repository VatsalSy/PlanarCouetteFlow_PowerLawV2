%% Getting the triple point contact line
% Author: Vatsal Sanjay
% vatsalsanjay@gmail.com
% Physics of Fluids
clc
clear
close
% 
% place = sprintf('lastNewt');
% ll=evalc(sprintf('!./getDatay %s', place));
% lolo=textscan(ll,'%f %f %f %f\n');
% yNewt=lolo{1}; uNewt=lolo{2}; shearNewt=lolo{3};
% yNewtTH = linspace(0.0,1.0,1001);
% uNewtTH = yNewtTH.*(1-yNewtTH/2);
% shearNewtTH = (1 - yNewtTH);
% 
% place = sprintf('lastShThn');
% ll=evalc(sprintf('!./getDatay %s', place));
% lolo=textscan(ll,'%f %f %f %f\n');
% yShTn=lolo{1}; uShTn=lolo{2}; shearShTn=lolo{3};
% yShTnTH = linspace(0.0,1.0,1001);
% uShTnTH = 0.5*(1/3 - (1-yShTnTH).^(3)/3);
% shearShTnTH = 0.5*(1-yShTnTH).^2;
% 
% place = sprintf('lastHB');
% ll=evalc(sprintf('!./getDatay %s', place));
% lolo=textscan(ll,'%f %f %f %f\n');
% yHB=lolo{1}; uHB=lolo{2}; shearHB=lolo{3};
% yHBTH = linspace(0.0,1.0,1001);
% tau_y = 0.25;
% Y = 1-tau_y;
% uHBTH = zeros(length(yHBTH),1);
% shearHBTH = zeros(length(yHBTH),1);
% for i = 1:1:length(yHBTH)
%     if (yHBTH(i) < Y)
%         y = yHBTH(i);
%         uHBTH(i) = (1/6)*(Y^3 - (Y-y).^3);
%         shearHBTH(i) = 0.5*(Y-y).^2;
%     else
%         uHBTH(i) = (1/6)*Y^3;
%         shearHBTH(i) = 0.0;
%     end
% end
% 
% 
% place = sprintf('lastBing');
% ll=evalc(sprintf('!./getDatay %s', place));
% lolo=textscan(ll,'%f %f %f %f\n');
% yBing=lolo{1}; uBing=lolo{2}; shearBing=lolo{3};
% yBingTH = linspace(0.0,1.0,1001);
% tau_y = 0.25;
% Y = 1-tau_y;
% uBingTH = zeros(length(yBingTH),1);
% shearBingTH = zeros(length(yBingTH),1);
% for i = 1:1:length(yBingTH)
%     if (yBingTH(i) < Y)
%         y = yBingTH(i);
%         uBingTH(i) = y.*(1-0.5*y) - tau_y*y;
%         shearBingTH(i) = (1-y) - tau_y;
%     else
%         uBingTH(i) = Y.*(1-0.5*Y) - tau_y*Y;
%         shearBingTH(i) = 0.0;
%     end
% end
% figure1 = figure('visible','on','WindowState','fullscreen','Color',[1 1 1]);
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% nskip = 4;
% plot(yNewtTH,uNewtTH,'k-','LineWidth',3,'DisplayName','(1.0, 0.0, 1.0): $u_{theory}$')
% plot(yNewt(2:nskip:end),uNewt(2:nskip:end),'k*','MarkerSize',15,'LineWidth',2,'DisplayName','(1.0, 0.0, 1.0): $u_{Basilisk}$')
% plot(yNewtTH,shearNewtTH,'k--','LineWidth',3,'DisplayName','(1.0, 0.0, 1.0): $\dot{\gamma}_{theory}$')
% plot(yNewt(2:nskip:end),shearNewt(2:nskip:end),'k.','MarkerSize',30,'LineWidth',2,'DisplayName','(1.0, 0.0, 1.0): $\dot{\gamma}_{Basilisk}$')
% 
% plot(yShTnTH,uShTnTH,'r-','LineWidth',3,'DisplayName','(1.0, 0.0, 0.5): $u_{theory}$')
% plot(yShTn(2:nskip:end),uShTn(2:nskip:end),'r*','MarkerSize',15,'LineWidth',2,'DisplayName','(1.0, 0.0, 0.5): $u_{Basilisk}$')
% plot(yShTnTH,shearShTnTH,'r--','LineWidth',3,'DisplayName','(1.0, 0.0, 0.5): $\dot{\gamma}_{theory}$')
% plot(yShTn(2:nskip:end),shearShTn(2:nskip:end),'r.','MarkerSize',30,'LineWidth',2,'DisplayName','(1.0, 0.0, 0.5): $\dot{\gamma}_{Basilisk}$')
% 
% plot(yHBTH,uHBTH,'g-','LineWidth',3,'DisplayName','(1.0, 0.25, 0.5): $u_{theory}$')
% plot(yHB(2:nskip:end),uHB(2:nskip:end),'g*','MarkerSize',15,'LineWidth',2,'DisplayName','(1.0, 0.25, 0.5): $u_{Basilisk}$')
% plot(yHBTH,shearHBTH,'g--','LineWidth',3,'DisplayName','(1.0, 0.25, 0.5): $\dot{\gamma}_{theory}$')
% plot(yHB(2:nskip:end),shearHB(2:nskip:end),'g.','MarkerSize',30,'LineWidth',2,'DisplayName','(1.0, 0.25, 0.5): $\dot{\gamma}_{Basilisk}$')
% 
% plot(yBingTH,uBingTH,'b-','LineWidth',3,'DisplayName','(1.0, 0.25, 1.0): $u_{theory}$')
% plot(yBing(2:nskip:end),uBing(2:nskip:end),'b*','MarkerSize',15,'LineWidth',2,'DisplayName','(1.0, 0.25, 1.0): $u_{Basilisk}$')
% plot(yBingTH,shearBingTH,'b--','LineWidth',3,'DisplayName','(1.0, 0.25, 1.0): $\dot{\gamma}_{theory}$')
% plot(yBing(2:nskip:end),shearBing(2:nskip:end),'b.','MarkerSize',30,'LineWidth',2,'DisplayName','(1.0, 0.25, 1.0): $\dot{\gamma}_{Basilisk}$')
% 
% box(axes1,'on');
% set(axes1,'DataAspectRatio',[1 1 1],'FontName','times new roman','FontSize',...
%     30,'FontWeight','bold','LineWidth',3,'PlotBoxAspectRatio',[1 1 1]);
% xlim([0.0 1.0])
% % ylim([0 1.0])
% axis square
% xlabel('\boldmath{$y$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
%             'FontName','times new roman',...
%             'Interpreter','latex');
% ylabel('\boldmath{$u, \dot{\gamma}$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
%     'FontName','times new roman',...
%     'Interpreter','latex'); 
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.767051496959868 0.341841841841842 0.209436598278227 0.588088088088088],...
%     'LineWidth',3,...
%     'Interpreter','latex',...
%     'FontSize',30,...
%     'EdgeColor',[1 1 1]);
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.786714285714287 0.911138274236123 0.107333333333332 0.0649895178197073],...
%     'String','\boldmath{$\left(\mu_0, \tau_y, n\right)$}',...
%     'Interpreter','latex',...
%     'FontSize',35,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.423619047619054 0.760797342192689 0.271619047619044 0.0735807267180653],...
%     'String','\boldmath{$\mu_{eq} = \frac{\tau_y}{2\|D_{ij}\|} + \mu_0\|D_{ij}\|^{n-1}$}',...
%     'Interpreter','latex',...
%     'FontSize',40,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.482547619047627 0.701701701701702 0.198999999999992 0.061605296137981],...
%     'String','\boldmath{$\tau = \tau_y + 2\mu_0D_{ij}^{n}$}',...
%     'Interpreter','latex',...
%     'FontSize',40,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.619452380952391 0.528528528528528 0.139476190476181 0.061605296137981],...
%     'String','\textbf{Newtonian}',...
%     'Interpreter','latex',...
%     'FontSize',40,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.642666666666677 0.347347347347347 0.116261904761894 0.061605296137981],...
%     'Color',[0 0 1],...
%     'String','\textbf{Bingham}',...
%     'Interpreter','latex',...
%     'FontSize',40,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.759928571428584 0.207207207207207 0.116261904761894 0.0616052961379812],...
%     'Color',[1 0 0],...
%     'String','\textbf{Shear Thinning}',...
%     'Interpreter','latex',...
%     'FontSize',40,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.759928571428586 0.138138138138138 0.116261904761894 0.0616052961379811],...
%     'Color',[0 1 0],...
%     'String','\textbf{Herschel-Bulkley}',...
%     'Interpreter','latex',...
%     'FontSize',40,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

%% For Bingham
place = sprintf('lastBing');
ll=evalc(sprintf('!./getDatay %s', place));
lolo=textscan(ll,'%f %f %f %f\n');
yBing=lolo{1}; uBing=lolo{2}; shearBing=lolo{3}; D2Bing=lolo{4};

yBingTH = linspace(0.0,1.0,1001);
tau_y = 0.25;
Y = 1-tau_y;
uBingTH = zeros(length(yBingTH),1);
shearBingTH = zeros(length(yBingTH),1);
for i = 1:1:length(yBingTH)
    if (yBingTH(i) < Y)
        y = yBingTH(i);
        uBingTH(i) = y.*(1-0.5*y) - tau_y*y;
        shearBingTH(i) = (1-y) - tau_y;
    else
        uBingTH(i) = Y.*(1-0.5*Y) - tau_y*Y;
        shearBingTH(i) = 0.0;
    end
end

figure1 = figure('visible','on','WindowState','fullscreen','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
DIJBing = D2Bing/sqrt(2);
DIJBingTH = 0.5*shearBingTH;
nskip = 4;
plot(yBing(2:nskip:end),uBing(2:nskip:end),'k*','MarkerSize',15,'LineWidth',2,'DisplayName','Basilisk: $u$')
plot(yBingTH,uBingTH,'k-','LineWidth',3,'DisplayName','Theory: $u$')

plot(yBing(2:nskip:end),shearBing(2:nskip:end),'r.','MarkerSize',30,'DisplayName','Basilisk: $\dot{\gamma}$')
plot(yBingTH,shearBingTH,'r-','LineWidth',3,'DisplayName','Theory: $\dot{\gamma}$')

plot(yBing(2:nskip:end),DIJBing(2:nskip:end),'b+','MarkerSize',15,'LineWidth',2,'LineWidth',3,'DisplayName','Basilisk: $\|D_{ij}\|$')
plot(yBingTH,DIJBingTH,'b-','LineWidth',3,'DisplayName','Theory: $\|D_{ij}\|$')

box(axes1,'on');
set(axes1,'DataAspectRatio',[1 1 1],'FontName','times new roman','FontSize',...
    30,'FontWeight','bold','LineWidth',3,'PlotBoxAspectRatio',[1 1 1]);
xlim([0.0 1.0])
axis square
xlabel('\boldmath{$y$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
            'FontName','times new roman',...
            'Interpreter','latex');
ylabel('\boldmath{$u, \dot{\gamma}, \|D_{ij}\|$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
    'FontName','times new roman',...
    'Interpreter','latex'); 
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.767646735055106 0.455955955955956 0.209436598278227 0.441941941941942],...
    'LineWidth',3,...
    'Interpreter','latex',...
    'FontSize',30,...
    'EdgeColor',[1 1 1]);
annotation(figure1,'textbox',...
    [0.481952380952388 0.762799344194691 0.203761904761898 0.0735807267180653],...
    'String','\boldmath{$\mu_{eq} = \frac{\tau_y}{2\|D_{ij}\|} + \mu_0$}',...
    'Interpreter','latex',...
    'FontSize',40,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.482547619047627 0.701701701701702 0.198999999999992 0.061605296137981],...
    'String','\boldmath{$\tau = \tau_y + 2\mu_0D_{ij}$}',...
    'Interpreter','latex',...
    'FontSize',40,...
    'FitBoxToText','off',...
    'EdgeColor','none');
