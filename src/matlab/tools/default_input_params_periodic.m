function params = default_input_params_periodic(name)
%DEFAULT_INPUT_PARAMS Creates a basic input structure that can be used to
%create domains. Fields are:
% name: simulation name
% periodic: flag for periodic domains
% panels: number of panels, can be single integer in which case all walls
% will be discretized with the same number of panels, or a vector of size
% #walls x 1
% plot_domain: flag to plot the domain along with normal vectors
% test: flag to test problem by applying special boundary conditions
% pressure_drop_x: pressure drop in x-direction
% pressure_drop_y: pressure drop in y-direction
% gmres_tol: GMRES tolerance
% eta: scaling factor in front of SLP in combined layer formulation

params = struct('name', name,...
                      'periodic', 1,...
                      'panels', 20,...
                      'plot_domain', 0,...
                      'test', 0,...
                      'pressure_drop_x', 1,...
                      'pressure_drop_y', 0,...
                      'box_size', [1,1],...
                      'gmres_tol', 1e-12,...
                      'eta', 1);


