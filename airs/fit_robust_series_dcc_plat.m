function fit_robust_series_dcc_plat(lati)
% FIT_ROBUST_SERIES_PLAT 
%
% 

addpath('/asl/packages/time')  % tai2dtime
addpath('/home/strow/Matlab/Math')  % remove_3sigma
addpath('/home/strow/Matlab/Breno/Scr/Math')  % Math_tsfit_...

% initialize output collection variables
all_b = zeros(1,2378,10);
all_berr = zeros(1,2378,10);
all_rms = zeros(1,2378);
anom = zeros(1,4509,2378);
all_bcorr = zeros(1,2378,10,10);
x1 = zeros(4509,1)*NaN;

% restore mission aggregate stats
% $$$ load('~/testoutput/save_test.mat')
load('~/testoutput/rates_temp/save_dcc.mat')

% convert rtime values to matlab datetime
%dmtime = tai2dtime(rtime_mean);
dmtime = tai2dnum(rtime_mean);  % can we make this work with datetime?

trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');

% load ig array from lls mat file
load('/home/strow/rates_12year.mat', 'ig');

warning('off','all');
    lati
% $$$     it = find(dmtime(:,lati) >= datetime('2003-05-30') & ...
% $$$               dmtime(:,lati) <= datetime('2015-05-30'));
    it = find(dmtime(:,lati) >= datenum('2003-05-30', 'yyyy-mm-dd') & ...
              dmtime(:,lati) <= datenum('2015-05-30', 'yyyy-mm-dd'));

    if numel(it) ~= 0
        for i=ig  % loop over channels (directly from array of channel values)
                  %x= seconds(dmtime(it,lati)-dmtime(it(1),lati));
            x = dmtime(it,lati)-dmtime(it(1),lati);
            %       y=squeeze(btobs(it,lati,i));
            %        y=squeeze(btclrcal(it,lati,i));
            %        y=squeeze(btcldcal(it,lati,i));
            %       y=squeeze(btobs(it,lati,i))-squeeze(btcldcal(it,lati,i));
% $$$         y=squeeze(btclrcal(it,lati,i))-squeeze(btobs(it,lati,i));
            y=squeeze(btcal(it,lati,i))-squeeze(btobs(it,lati,i));
            ii = remove_3sigma(y);
            x = x(ii);y=y(ii);

            % may need to wrap this in try-catch and set variables as 0
            % or NaN in event of failure 
            [b stats] = Math_tsfit_lin_robust(x,y,4);
            %      [b stats] =
            %      Math_tsfit_lin_robust(dmtime(it,lati)-dmtime(it(1),lati),squeeze(btcldcal(it,lati,i)),4);
            all_b(1,i,:) = b;
            all_rms(1,i) = stats.s;
            x1=x1*NaN;
            x1(ii) = stats.resid;
            anom(1,:,i) = x1;
            all_berr(1,i,:) = stats.se;
            all_bcorr(1,i,:,:) = stats.coeffcorr;
        end  % end for
    end  % end if
warning('on','all');
rates_clr = squeeze(all_b(:,:,2));
rates_clr_err = squeeze(all_berr(:,:,2));
outfile = fullfile('/home/sbuczko1/testoutput/rates_temp', ['rates_' int2str(lati) '_dcc.mat']);
save(outfile, 'rates_clr', 'rates_clr_err', 'ig')


%% ****end function fit_robust_series_plat****