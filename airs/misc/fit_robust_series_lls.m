addpath('/asl/packages/time')  % tai2dtime
addpath('/home/strow/Matlab/Math')  % remove_3sigma
addpath('/home/strow/Matlab/Breno/Scr/Math')  % Math_tsfit_...

% initialize output collection variables
all_b = zeros(40,2378,10);
all_berr = zeros(40,2378,10);
all_rms = zeros(40,2378);
anom = zeros(40,4509,2378);
all_bcorr = zeros(40,2378,10,10);

% restore mission aggregate stats
load('~/testoutput/save_dcc.mat')

% convert rtime values to matlab datetime
%dmtime = tai2dtime(rtime_mean);
dmtime = tai2dnum(rtime_mean);  % can we make this work with datetime?

trace.RunDate = datetime('now','TimeZone','local','Format', ...
                         'd-MMM-y HH:mm:ss Z');

% load ig array from lls mat file
load('/home/strow/rates_12year.mat', 'ig');

warning('off','all');
for lati = 1:40
    lati
    % select for May 30, 2003 to May 30, 2015 (could probably build
    % this once outside loop??)
% $$$     it = find(dmtime(:,lati) >= datetime('2003-05-30') & ...
% $$$               dmtime(:,lati) <= datetime('2015-05-30'));
    it = find(dmtime(:,lati) >= datenum('2003-05-30', 'yyyy-mm-dd') & ...
              dmtime(:,lati) <= datenum('2015-05-30', 'yyyy-mm-dd'));

    % dcc will have empty lat bins. check for this and move on
    if numel(it) == 0
        continue
    end
    
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
        all_b(lati,i,:) = b;
        all_rms(lati,i) = stats.s;
        x1 = zeros(4509,1)*NaN;
        x1(ii) = stats.resid;
        anom(lati,:,i) = x1;
        all_berr(lati,i,:) = stats.se;
        all_bcorr(lati,i,:,:) = stats.coeffcorr;
    end
end
warning('on','all');
rates_clr = squeeze(all_b(:,:,2));
rates_clr_err = squeeze(all_berr(:,:,2));
save ~/testoutput/2015/airs/rates/rates_dcc rates_clr rates_clr_err ig
