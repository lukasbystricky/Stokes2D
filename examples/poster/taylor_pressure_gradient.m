% Poster figures
close all
clearvars
clc

% create input structure
input_params = default_input_params('couette', 0);

% modify structure as needed, or add additional problem-dependent params
input_params.panels = [22,11];%40 20
input_params.plot_domain = 0;

% set up radii and centers, centers given as complex numbers (x+iy)
input_params.radii = [2, 1]; % the outer wall must be the first radii
input_params.centers = [0, 0];
input_params.omega = 1; % angular velocity of inner cylinder
input_params.slip = 0;

problem = couette(input_params);

% solve the problem
z = problem.domain.z;

rhs = [real(problem.boundary_conditions(z));... 
       imag(problem.boundary_conditions(z))];
solution = solve_stokes(problem,rhs,'fmm',0);
solution.local_indices = 1:length(solution.q);

%% Off-surface
% grid in polar coordinates with M number of radial and angular points
M = 200;
r = linspace(input_params.radii(2)+1e-6, input_params.radii(1) - 1e-6, M);
h = (2*pi)/M;
theta = (0:h:2*pi);
[R, T] = meshgrid(r, theta);

X = R.*cos(T);
Y = R.*sin(T);

[Pxc, Pyc, Px, Py, ~, ~] = evaluate_pressure_gradient(solution, X, Y);

%% Plot
close all;
save_fig = 0;
pub_fig = 0;

c1 = [161 63 54];
c2 = [212 123 93];
c3 = [214 176 155];
c4 = [170 179 194];
c5 = [91 118 161];
c6 = [29 51 101];

R = [29 91 170 214 212 161]./255*1.15;
G = [51 118 179 176 123 63]./255*1.15;
B = [101 161 194 155 93 54]./255*1.15;
M = length(R);

E_reg_quad = abs(Px+1i*Py);
E_sq = abs(Pxc+1i*Pyc);

E_reg_quad(E_reg_quad < 1e-15) = 1e-15;
E_sq(E_sq < 1e-15) = 1e-15;

sfigure(1);
if pub_fig
    publication_fig;
end
plot(problem.domain.z,'k.');
hold on
[C,h]=contourf(X,Y,log10(E_reg_quad),M);
set(h,'LineColor','none');
xlabel(colorbar(), 'log_{10} E_{abs}');
caxis([-15 0]);
colormap([R(:), G(:), B(:)]);
%ax = gca;
%ax.XLim = [0 2];
%ax.YLim = [0 2];
axis off;
axis equal;

sfigure(2);
if pub_fig
    publication_fig;
end
plot(problem.domain.z,'k.');
hold on
[c2,h2]=contourf(X,Y,log10(E_sq),M);
set(h2,'LineColor','none');
xlabel(colorbar(), 'log_{10} E_{abs}');
caxis([-15 0]);
colormap([R(:), G(:), B(:)]);
axis off;
axis equal;
%ax = gca;
%ax.XLim = [0 2];
%ax.YLim = [0 2];

if save_fig
    printpdf(1,'taylor_couette_reg_quad');
    printpdf(2,'taylor_couette_sq');
end

function publication_fig
set(gcf,'PaperPositionMode','auto')
P = get(gcf,'Position');
set(gcf,'Position',[P(1) P(2) 275 220])
set(gca,'FontName','Times','FontSize',8)
box off
end

function printpng(num, name)
    sfigure(num);
    fname = [name, '.png'];
    print('-dpng', '-r200', fname)
    disp(['Wrote ' fname])
    system(['mogrify -trim ' fname]);
end    

function printpdf(num, name)
    sfigure(num);
    fname = [name, '.pdf'];
    exportgraphics(gcf,fname,'ContentType','vector')
    disp(['Wrote ' name])
    %system(['mogrify -trim ' fname]);
end

function h = sfigure(f)
if nargin==0
    h = figure();
elseif ishandle(f)
    set(0, 'CurrentFigure', f);
else
    h = figure(f);
end
end
