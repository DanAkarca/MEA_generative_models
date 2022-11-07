%% supplementary figure 3
% written by danyal akarca
clear; clc;
% set directory
gendir = '/imaging/astle/users/da04/PhD/hd_gnm_generative_models/raw_data/ToShare_Aug2022/';
% 50k
load(strcat(gendir,'50k_rodent_dev/M03212_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat')); div7 = gm_data;
load(strcat(gendir,'50k_rodent_dev/M03212_rec2_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div10 = gm_data;
load(strcat(gendir,'50k_rodent_dev/M03212_rec3_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div12 = gm_data;
load(strcat(gendir,'50k_rodent_dev/M03212_rec4_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div14 = gm_data;
rodent_50 = {div7 div10 div12 div14};
% 100k
load(strcat(gendir,'100k_rodent_dev/100k_combined_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat')); div14 = gm_data;
load(strcat(gendir,'100k_rodent_dev/100k_combined_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div21 = gm_data;
load(strcat(gendir,'100k_rodent_dev/100k_combined_rec1_min_rate_0.01hz_alpha_0.001_lag10_jitter10_prob.mat'));  div28 = gm_data;
rodent_100 = {div14 div21 div28}; 
% human iPSCs
load(strcat(gendir,'human_ipsc_gn/GN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat')); gn = gm_data;
load(strcat(gendir,'human_ipsc_mn/MN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat'));  mn = gm_data;
load(strcat(gendir,'human_ipsc_dn/DN_gm_data_min_r_0.01_alpha_0.001_1800s_min_rate_001hz_lag10_jitter10_prob.mat'));  dn = gm_data;
ipsc = {gn mn dn};
% human COs
load(strcat(gendir,'human_organoids/M03912_GNM_data_min_rate_0.01hz_0.001_alpha_1800s_min_rate_001hz_lag10_jitter10_prob.mat')); hco = gm_data;
%% plot the euclidean distances
% set the scale
maxscale = 4300;
% initialise
m = {};

% 50k
% settings
div = 4;
nculture = 6;
% matrices
h = figure; h.Position = [100 100 1200 155];
for culture = 1:nculture;
    subplot(1,nculture,culture);
    D = squareform(pdist(rodent_50{div}.table_xy{culture}));
    imagesc(D);
    b = gca; b.TickDir = 'out'; axis off;
    % keep the maximum to then scale
    m{1}(culture) = max(D,[],'all');
    % keep the upper triangle
    ind = find(triu(D,1));
    cost{1}{culture} = D(ind);
    % limits
    caxis([0 maxscale]);
end
% histogram
h = figure; h.Position = [100 100 650 350];
for culture = 1:nculture;
    [y x] = ksdensity(cost{1}{culture});
    plot(x,y,...
        'linewidth',4);
    hold on;
end
xlabel('Euclidean (\mum)'); ylabel('Frequency'); yticks([]);
b = gca; b.TickDir = 'out'; xlim([0 inf]);
box off; b.FontName = 'Arial'; b.FontSize = 16;
sgtitle('n=6 Rodent 50k cultures (DIV14)','FontSize',16);
xlim([0 maxscale]);

% 100k
% settings
div = 3;
nculture = 12;
% matrices
h = figure; h.Position = [100 100 1200 350];
for culture = 1:nculture;
    subplot(2,nculture/2,culture);
    D = squareform(pdist(rodent_100{div}.table_xy{culture}));
    imagesc(D);
    b = gca; b.TickDir = 'out'; axis off;
    % keep the maximum to then scale
    m{2}(culture) = max(D,[],'all');
    % keep the upper triangle
    ind = find(triu(D,1));
    cost{2}{culture} = D(ind);
    % limits
    caxis([0 maxscale]);
end
% histogram
h = figure; h.Position = [100 100 650 350];
for culture = 1:nculture;
    [y x] = ksdensity(cost{2}{culture});
    plot(x,y,...
        'linewidth',4);
    hold on;
end
xlabel('Euclidean (\mum)'); ylabel('Frequency'); yticks([]);
b = gca; b.TickDir = 'out'; xlim([0 inf]);
box off; b.FontName = 'Arial'; b.FontSize = 16;
sgtitle('n=12 Rodent 100k cultures (DIV28)','FontSize',16);
xlim([0 maxscale]);

% gn
% settings
nculture = 8;
% matrices
h = figure; h.Position = [100 100 1200 350];
for culture = 1:nculture;
    subplot(2,6,culture);
    D = squareform(pdist(gn.table_xy{culture}));
    imagesc(D);
    b = gca; b.TickDir = 'out'; axis off;
    % keep the maximum to then scale
    m{3}(culture) = max(D,[],'all');
    % keep the upper triangle
    ind = find(triu(D,1));
    cost{3}{culture} = D(ind);
    % limits
    caxis([0 maxscale]);
end
% histogram
h = figure; h.Position = [100 100 650 350];
for culture = 1:nculture;
    [y x] = ksdensity(cost{3}{culture});
    plot(x,y,...
        'linewidth',4);
    hold on;
end
xlabel('Euclidean (\mum)'); ylabel('Frequency'); yticks([]);
b = gca; b.TickDir = 'out'; xlim([0 inf]);
box off; b.FontName = 'Arial'; b.FontSize = 16;
sgtitle('n=8 Human 100k glutamatergic iPSC cultures (DIV28)','FontSize',16);
xlim([0 maxscale]);

% mn
% settings
nculture = 7;
% matrices
h = figure; h.Position = [100 100 1200 350];
for culture = 1:nculture;
    subplot(2,6,culture);
    D = squareform(pdist(mn.table_xy{culture}));
    imagesc(D);
    b = gca; b.TickDir = 'out'; axis off;
    % keep the maximum to then scale
    m{4}(culture) = max(D,[],'all');
    % keep the upper triangle
    ind = find(triu(D,1));
    cost{4}{culture} = D(ind);
    % limits
    caxis([0 maxscale]);
end
% histogram
h = figure; h.Position = [100 100 650 350];
for culture = 1:nculture;
    [y x] = ksdensity(cost{4}{culture});
    plot(x,y,...
        'linewidth',4);
    hold on;
end
xlabel('Euclidean (\mum)'); ylabel('Frequency'); yticks([]);
b = gca; b.TickDir = 'out'; xlim([0 inf]);
box off; b.FontName = 'Arial'; b.FontSize = 16;
sgtitle('n=7 Human 100k motor iPSC cultures (DIV28)','FontSize',16);
xlim([0 maxscale]);

% dn
% settings
nculture = 6;
% matrices
h = figure; h.Position = [100 100 1200 155];
for culture = 1:nculture;
    subplot(1,6,culture);
    D = squareform(pdist(dn.table_xy{culture}));
    imagesc(D);
    b = gca; b.TickDir = 'out'; axis off;
    % keep the maximum to then scale
    m{5}(culture) = max(D,[],'all');
    % keep the upper triangle
    ind = find(triu(D,1));
    cost{5}{culture} = D(ind);
    % limits
    caxis([0 maxscale]);
end
% histogram
h = figure; h.Position = [100 100 650 350];
for culture = 1:nculture;
    [y x] = ksdensity(cost{5}{culture});
    plot(x,y,...
        'linewidth',4);
    hold on;
end
xlabel('Euclidean (\mum)'); ylabel('Frequency'); yticks([]);
b = gca; b.TickDir = 'out'; xlim([0 inf]);
box off; b.FontName = 'Arial'; b.FontSize = 16;
sgtitle('n=6 Human 100k dopaminergic iPSC cultures (DIV28)','FontSize',16);
xlim([0 maxscale]);

% hco
% settings
nculture = 6;
% matrices
h = figure; h.Position = [100 100 1200 155];
for culture = 1:nculture;
    subplot(1,6,culture);
    D = squareform(pdist(hco.table_xy{culture}));
    imagesc(D);
    b = gca; b.TickDir = 'out'; axis off;
    % keep the maximum to then scale
    m{6}(culture) = max(D,[],'all');
    % keep the upper triangle
    ind = find(triu(D,1));
    cost{6}{culture} = D(ind);
    % limits
    caxis([0 maxscale]);
end
% histogram
h = figure; h.Position = [100 100 650 350];
for culture = 1:nculture;
    [y x] = ksdensity(cost{6}{culture});
    plot(x,y,...
        'linewidth',4);
    hold on;
end
xlabel('Euclidean (\mum)'); ylabel('Frequency'); yticks([]);
b = gca; b.TickDir = 'out'; xlim([0 inf]);
box off; b.FontName = 'Arial'; b.FontSize = 16;
sgtitle('n=6 Human cerebral organoids','FontSize',16);
xlim([0 maxscale]);