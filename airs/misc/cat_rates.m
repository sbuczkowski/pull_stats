load(fullfile('~/testoutput/rates_temp', 'rates_1.mat'))

rates_clr_temp = rates_clr;
rates_clr_err_temp = rates_clr_err;
clear rates_clr rates_clr_err

for i=2:40
    load(fullfile('~/testoutput/rates_temp', ['rates_' int2str(i) '.mat']))
    rates_clr_temp = cat(1, rates_clr_temp, rates_clr);
    rates_clr_err_temp = cat(1, rates_clr_err_temp, rates_clr_err);
    clear rates_clr rates_clr_err
end

rates_clr = rates_clr_temp;
rates_clr_err = rates_clr_err_temp;

clear i rates_clr_temp rates_clr_err_temp

save(fullfile('~/testoutput/2015/airs/rates', 'rates_rand_btobs.mat'), 'rates_clr', 'rates_clr_err', 'ig')
    