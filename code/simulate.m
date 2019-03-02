close all;

%% DECLARE PRE-SET GLOBAL PARAMETERS
global F_g_set; % gravity (1 x 3)
global k_v;
global k_of;
global k_r;
global k_a;
global k_air;
global m; % number of bubbles
global R; % radii (m x 1)
global PO;
global RO;
global mo;
global k_ro;
global k_ao;

%% SET PARAMETERS
% <PARAMETERS OF ENVIRONMENT> =============================================
F_g_set = [0,-10,-10]; % gravity (1 x 3)
k_air = 5;% friction coefficient of air
ax = 100;% set axis limits

% <PARAMETERS OF OBJECTS> =================================================
k_ro = 5*1000; % repulse coefficient of objects
k_ao = 50; % attract coefficient of objects
k_of = 5; % surface friction coefficient of objects

% NOTE THAT THE FLOOR IS SET AS A LARGE SPHERE OBJECTS

% % MODE_1: no objects(/obstacles)
% PO = [0,0,-10000;0,-10100,0];
% RO = [10000;10000];
% % MODE_2: cylinder object
% PO = [0,0,-10000;0,-10100,0;0,0,0;0,0,2;0,0,4;0,0,6;0,0,8];
% RO = [10000;10000;20;20;20;20;20];
% MODE_3: hemisphere object
PO = [0,0,-10000;0,-10100,0;0,0,0];
RO = [10000;10000;20];

mo = length(RO);% number of objects

% <PARAMETERS OF BUBBLES> =================================================
k_v = 50;% viscosity coefficient of liquid
k_r = 5; % repulse coefficient of bubbles
k_a = 5; % attract coefficient of bubbles
m = 100; % number of bubbles

% MODE_1: set generate rigion: from one side
xa = -5;
xb = 5;
ya = 80;
yb = 100;
za = 20;
zb = 40;
% % MODE_2: set generate rigion: from above
% xa = -5;
% xb = 5;
% ya = -5;
% yb = 5;
% za = 40;
% zb = 60;

% inputs
X = xa + (xb-xa).*rand(m,1);
Y = ya + (yb-ya).*rand(m,1);
Z = za + (zb-za).*rand(m,1);
P = [X,Y,Z]; % positons (m x 3)
V = zeros(m,3); % velocities (m x 3)
R = 5.* ones(m,1); % radii (m x 1)

%% START NUMERICAL ITERATION

for i = 1:300
    [P,V] = computeNext(P, V); % compute next Velocity & Position of bubbles
%     plot3(P(:,1), P(:,2), P(:,3) ,'o','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerFaceColor','r')
    scatter3(P(:,1),P(:,2),P(:,3),100,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.1)
    
    hold on
    
    % draw semisphere
    phi=0:0.01:pi;
    Theta=0:0.02:2.1*pi;
    x=RO(3)*sin(phi)'*cos(Theta);
    y=RO(3)*sin(phi)'*sin(Theta);
    z=(RO(3)^2-x.^2-y.^2).^(1/2);
    surf(x,y,z,'edgecolor','none','FaceAlpha',0.4)
    
%     % draw cylinder
%     [x,y,z]=cylinder(RO(3),100);
%     z=ax*z;
%     surf(x,y,z,'edgecolor','none','FaceAlpha',0.4)

    % draw the floor
    x=-ax:50:ax;
    y=x;
    [x,y]=meshgrid(x,y);
    z=x*0;
    surf(x,y,z,'edgecolor','none','FaceAlpha',0.6)
    colormap([0 0 0])
    
    axis([-ax,ax,-ax,ax,-ax,ax])
    grid on
    hold off
    drawnow
end
