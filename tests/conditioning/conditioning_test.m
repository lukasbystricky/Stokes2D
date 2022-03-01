close all
clearvars
clc

% create input structure
input_params = default_input_params('pouseuille_demo', 1);

% modify structure as needed
input_params.box_size = [2,5];
input_params.h = 0.5;    % pipe walls at +-0.5
input_params.panels = 5;
input_params.eta = 1;
input_params.plot_domain = 0;
input_params.slip = 1;
input_params.A = 1e-10;
input_params.d = 0.5;

alphavec = [0, 1e-8];
M = length(alphavec);
input_params_cell = cell(M,1);

for i = 1:M
    input_params.alpha = alphavec(i);
    input_params_cell{i} = input_params;
end

%%
rankvec = zeros(M,1);
condvec = zeros(M,1);
K = zeros(2*input_params.panels*16*2+2,2*input_params.panels*16*2+2,M);

for i = 1:M
    disp([i M]);
    problem = flat_pipe_periodic(input_params_cell{i});
    %problem = sine_pipe_periodic(input_params_cell{i});
    N = length(problem.domain.z);
    A = zeros(2*N+2);
    X = eye(2*N+2);
    parfor j = 1:2*N+2
        A(:,j) = matvec_combined_slip(X(:,j),problem);
    end
    rankvec(i) = rank(A);
    condvec(i) = cond(A);
    K(:,:,i) = A;
    
    %figure;
    %plot(eig(A),'x');
    %grid on;
    %title(num2str(alphavec(i)));
end
rankvec
condvec

%%
diff = K(1:end-2,1:end-2,1)-K(1:end-2,1:end-2,2);
figure;
mesh(log10(abs(diff)));
title('log10(abs(diff))');

% %% 
% for i = 1:M
%     A = K(:,:,i);
%     figure;
%     mesh(A);
%     figure;
%     spy(A);
% end

%%
figure;
loglog(alphavec,condvec,'x--');
xlabel('alpha');
ylabel('condition number');
title('top sine bottom flat wall, eta=1');
grid on;

% %% plot velocity slip profiles
% p = -0.5;
% h = 0.5;
% x = linspace(0,0.5,100);
% y = linspace(-0.5,0.5,100);
% alphavec = [0; 0.25; 0.5; 0.8; 1];
% 
% exact_solution_u = @(y,alpha) p/2*(y.^2-h^2) - h*p*alpha;
% 
% figure;
% hold on;
% for i = 1:length(alphavec)
%     plot(exact_solution_u(y,alphavec(i)),y);
% end
% xlim([0,0.5]);
% ylim([-0.5,0.5]);
% title('Pouseille flow velocity profile for slip');
% legend('alpha=0','alpha=0.25','alpha=0.5','alpha=0.8','alpha=1.0');
% xlabel('Horizontal velocity u');
% ylabel('Vertical position y');
% grid on;
