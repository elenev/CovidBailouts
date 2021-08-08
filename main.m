% Master file

% Run this file to reproduce all results in the paper

% Results folder will be filled with
% - PDF and EPS versions of all graphs
% - simres.mat file containing a hashmap storing every 
%   statistic cited in tables and text of the paper

%% Settings
no_par_processes = feature('numcores');
batchmode = true;

%% Solve and simulate model
disp('Solving models.')
expers = {'xi88midxigrid', 'xi88pandemic'};

for exper = expers 
   % Initial round
   disp('Economy: 
   expername = [exper{:},'_ini0'];
   maxit = 100;
   guess_mode = 'no_guess';
   outfname=['env_',expername,'.mat'];
   
   main_create_env;
   
   exper_path = outfname;
   outname = ['res_',expername,'.mat'];
   price_zns = false;
   
   main_run_exper;
   
   % Fine grid round
   expername = exper{:};
   maxit = 30;
   guess_mode = 'guess';
   guess_path = outname;
   outfname=['env_',expername,'.mat'];
   
   main_create_env;
   
   exper_path = ['env_',expername,'.mat'];
   outname = ['res_',expername,'.mat'];
   price_zns = true;
   
   main_run_exper;
   
   % Simulate
   resfile = ['res_',expername];
   
   sim_stationary;
   
   % Anticipated (+ prob ) shock IRFs
   sim_trans_cluster;
end

%% Simulate MIT shocks
% old -> new ("New Normal")
resfile = 'res_xi88pandemic';
start_resfile = 'res_xi88midxigrid';
sim_trans_mit;

% old -> old (Appendix)
resfile = 'res_xi88midxigrid';
start_resfile = 'res_xi88midxigrid';
sim_trans_mit;

%% Comparative welfare measure
welfare(10);

%% Make plots
compare_MIT_shocks_2;
compare_MIT_shocks_components;

%% Store stats
writeStats;