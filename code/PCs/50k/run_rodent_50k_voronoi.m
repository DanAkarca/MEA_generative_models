%% Run voronoi tesselation of 50k rodent primary cortical neuronal cultures
% Voronoi tesselation code written by Betzel Lab: https://www.brainnetworkslab.com/coderesources
% Adapted for use by Dr Danyal Akarca, University of Cambridge, 2022
%% load data for rodent primary cortical neurons
function output = run_rodent_50k_voronoi(pipeline,network,div);
%% add relevant paths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
addpath('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/voronoi');
%% load sttc data and set save directory
% change to data directory
cd('/imaging/astle/users/da04/PhD/hd_gnm_generative_models/ForAlex_updated/Oct_2021/Euler_50k_rodent/tracking/');
% load sttc data
a = load('gm_data_shared_div7_10_12_14_min_r_0.1_alpha_0.001_lag5_jitter5.mat');
b = load('gm_data_shared_div7_10_12_14_min_r_0.1_alpha_0.001_lag10_jitter10.mat');
c = load('gm_data_shared_div7_10_12_14_min_r_0.1_alpha_0.001_lag20_jitter20.mat');
% collect data
data = {a.gm_data b.gm_data c.gm_data};
% set save directory
sdir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/rodent_50k_qc/';
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
Atgt = data{pipeline}.table_bu_shared{network}{div};
% take the euclidean
D = squareform(pdist(data{pipeline}.table_xy_shared{network}{div}));
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
    disp(sprintf('running network_%g_%g_%g model %s...',pipeline,network,div,modeltypes(model)));
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
save(sprintf('rodent_50k_qc_%g_%g_%g_generative_model.mat',pipeline,network,div),'output','-v7.3');
end
