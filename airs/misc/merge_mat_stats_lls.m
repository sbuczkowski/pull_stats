%function  = merge_mat_stats_lls(statpath, regexp)
% MERGE_MAT_STATS_LLS 
%
% 

a = dir(fullfile(statpath,regexp));
load(fullfile(statpath, a(1).name));

for i=2:length(a)
   g = load(fullfile(statpath, a(i).name));
   %    bias = [bias; g.bias];
   % bias_std_clr = [bias_std_clr; g.bias_std_clr];
   % btcldcal = [btcldcal; g.btcldcal];
   % btclrcal = [btclrcal; g.btclrcal];
   btcal = [btcal; g.btcal];
   btobs = [btobs; g.btobs];
   %    btcal = [btcal; g.btcal];
   %    btobs = [btobs; g.btobs];
   %    bias_std = [bias_std; g.bias_std];
   lat_mean = [lat_mean; g.lat_mean];
   lon_mean = [lon_mean; g.lon_mean];
   solzen_mean = [solzen_mean; g.solzen_mean];
   satzen_mean = [satzen_mean; g.satzen_mean];
   rtime_mean = [rtime_mean; g.rtime_mean];
   stemp_mean = [stemp_mean; g.stemp_mean];
   count = [count; g.count];
   %   iudef4_mean = [iudef4_mean g.iudef4_mean];
end


%% ****end function merge_mat_stats_lls****
