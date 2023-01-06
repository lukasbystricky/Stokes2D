close all;
clearvars;
clc;

load periodic_corrugated_channel_ellipses;
%rng(123445562); % was used

problem = solution.problem;
npan = 40;
nq = 16;
%vol = 2 - 0.393816218162103;  % volume (area) of domain

% off-surface velocity
[U,V,X,Y] = evaluate_velocity(solution, 200, 'fmm', 0, 'verbose', 0);

% on-surface vorticity
Omega_surf = evaluate_vorticity_on_surface(solution,solution,'fluid');

%% Plot
clear gca;
close all;
save_fig = 0;
pub_fig = 0;

R = [29 91 170 214 212 161]./255*1.15;
G = [51 118 179 176 123 63]./255*1.15;
B = [101 161 194 155 93 54]./255*1.15;

val = linspace(0,1,200);

sfigure(1);
if pub_fig
    publication_fig;
end
hold on;
for i = 1:22
    plot(problem.domain.z((i-1)*npan*nq+1:i*npan*nq),'k');
end
surfc(X,Y,abs(U+1i*V));
shading interp;
xlabel(colorbar(), '$\|\mathbf{u}\|$','interpreter','latex');
cmap = colormap([R(:), G(:), B(:)]);
hsv = rgb2hsv(cmap);
cm_data=interp1(linspace(0,1,size(cmap,1)),hsv,val);
cm_data=hsv2rgb(cm_data);
colormap(cm_data);
axis equal;
set(gca,'color','none');    % remove white background
set(gca,'XColor','none');   % remove x-axis completely
yticks([-1/2,1/2]);         % where to place y-ticks
yticklabels({'0','1'});     % labels at specific y-ticks

% Plot vorticity
% For the Stokes equations, on boundaries that have a no-slip boundary 
% condition, the shear stress is equivalent to the vorticity
input_params.n_periods_top = 4;
input_params.n_periods_bottom = 2;
input_params.amplitude_top = 0.05;
k1 = 2*pi*input_params.n_periods_top/problem.Lx;
k2 = 2*pi*input_params.n_periods_bottom/problem.Lx;
a1 = input_params.amplitude_top;
a2 = input_params.amplitude_top;

bottom_wall = @(T,h) geometry_periodic_channel(@(t) a2*sin(k2*t) - h, ...
            @(t) a2*k2*cos(k2*t), @(t) -a2*k2^2*sin(k2*t),T, -1, problem.Lx);

T = linspace(0,2*pi,1000);
%d = logspace(-5,-1,6);
d = [1e-2 1e-3 1e-4 1e-5];
hvec = 0.5 - d;

sfigure(2);
if pub_fig
    publication_fig;
end
hold on;
for i = 1:length(hvec)
    h = hvec(i);
    [Zline,Zpline,Zppline] = bottom_wall(T,h);
    Xline = real(Zline);
    Yline = imag(Zline);

    omega = evaluate_vorticity(solution, Xline, Yline, 'fmm', 0, 'verbose', 0);
    plot(T,omega,'color',[R(i) G(i) B(i)],'linewidth',0.5);
end

TT = linspace(0,2*pi,npan*nq);
plot(TT,Omega_surf(end-npan*nq+1:end),'color',[R(end) G(end) B(end)],'linewidth',0.5);
xlim([0,2*pi]);
xlabel('t');
ylabel('Vorticity');
legend('$d = 10^{-2}$','$d = 10^{-3}$','$d = 10^{-4}$','$d = 10^{-5}$','on-surface','location','southwest','interpreter','latex');

if save_fig
    saveplot(1,'periodic_corrugated_channel_velocity_abs');
    %printpdf(1,'periodic_corrugated_channel_velocity_abs');
    %printpdf(2,'periodic_corrugated_channel_vorticity_bottom_wall_near_on_surface');
end

%% Functions
function publication_fig
set(gcf,'PaperPositionMode','auto')
P = get(gcf,'Position');
set(gcf,'Position',[P(1) P(2) 275 220])
set(gca,'FontName','Times','FontSize',8)
box off
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

function printpdf(num, name)
    sfigure(num);
    fname = [name, '.pdf'];
    exportgraphics(gcf,fname,'ContentType','vector')
    disp(['Wrote ' name])
end

function saveplot(num, filename)
    sfigure(num);
    publication_fig();
    path = filename;
    print('-dpng', '-r1000', path)
    system(['mogrify -trim ' path])    
    disp(['Wrote ' path])
end
