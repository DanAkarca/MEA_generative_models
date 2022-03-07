%% Run voronoi tesselation of human iPSC neuronal cultures
% Voronoi tesselation code written by Betzel Lab: https://www.brainnetworkslab.com/coderesources
% Adapted for use by Dr Danyal Akarca, University of Cambridge, 2022
function output = run_human_100k_voronoi(pipeline,network,div,type);
%% add relevant paths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/voronoi');
%% load sttc data and set save directory
% change to data directory
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ForAlex_updated/Euler_100k_IPSC');
% load sttc data
g1 = load('GN_gm_data_min_r_0.1_alpha_0.001_lag5_jitter5_prob.mat');
g2 = load('GN_gm_data_min_r_0.1_alpha_0.001_lag10_jitter10_prob.mat');
g3 = load('GN_gm_data_min_r_0.1_alpha_0.001_lag20_jitter20_prob.mat');
m1 = load('MN_gm_data_min_r_0.1_alpha_0.001_lag5_jitter5_prob.mat');
m2 = load('MN_gm_data_min_r_0.1_alpha_0.001_lag10_jitter10_prob.mat');
m3 = load('MN_gm_data_min_r_0.1_alpha_0.001_lag20_jitter20_prob.mat');
d1 = load('DN_gm_data_min_r_0.1_alpha_0.001_lag5_jitter5_prob.mat');
d2 = load('DN_gm_data_min_r_0.1_alpha_0.001_lag10_jitter10_prob.mat');
d3 = load('DN_gm_data_min_r_0.1_alpha_0.001_lag20_jitter20_prob.mat');
% collect the data together
gluta = {g1.gm_data g2.gm_data g3.gm_data};
motor = {m1.gm_data m2.gm_data m3.gm_data};
dopa = {d1.gm_data d2.gm_data d3.gm_data};
data = {gluta motor dopa};
% set save directory
sdir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/human_100k_qc/';
%% initialise model information
% set number of models
nmodels     = 13;
% set model type
modeltypes = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','deg-prod'});
% set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'},{'powerlaw'}];
%% set up parameters for the tesselation
% set eta and gamma limits
eta = [-10,10];
gam = [-10,10];
% parameters related to the optimization
pow = 2; % severity
nlvls = 5; % number of steps
nreps = 4000; % number of repetitions/samples per step
%% run the tesselation procedure
% set target network to this pipeline
Atgt = data{type}{pipeline}.table_bu{network}{div};
% take the euclidean
D = squareform(pdist(data{type}{pipeline}.table_xy{network}{div}));
% take the nnode
nnode = size(D,1);
% get key observed statistics
x = cell(4,1);
x{1} = sum(Atgt,2);
x{2} = clustering_coef_bu(Atgt);
x{3} = betweenness_bin(Atgt)';
x{4} = D(triu(Atgt,1) > 0);
% set seed
Aseed = zeros(nnode,nnode);
% set number of connections
m = nnz(Atgt)/2;
% initialise
output = struct;
output.energy = zeros(nmodels,nlvls.*nreps);
output.ks = zeros(nmodels,nlvls.*nreps,4);
output.networks = zeros(nmodels,m,nlvls.*nreps);
output.parameters = zeros(nmodels,nlvls.*nreps,2);
% nnode
n = size(Atgt,1);
for model = 1:nmodels;
    % print text
    disp(sprintf('running network_%g_%g_%g_%g model %s...',pipeline,network,div,type,modeltypes(model)));
    % run the model
    [E,K,N,P] = fcn_sample_networks(Atgt,Aseed,D,m,modeltypes(model),modelvar,nreps,nlvls,eta,gam,pow);
    % store the output
    output.energy(model,:) = E;
    output.ks(model,:,:) = K;
    output.networks(model,:,:) = N;
    output.parameters(model,:,:) = P;
end
%% save the subject
% change directory
cd(sdir);
% save file
save(sprintf('human_100k_qc_%g_%g_%g_%g_generative_model.mat',pipeline,network,div,type),'output','-v7.3');
end