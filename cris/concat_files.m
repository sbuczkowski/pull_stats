rtp_addpaths

maxobs = 100000;
% $$$ maxobs = 65000;

days = dir('044*');
ndays = length(days);

fprintf(1, '%d days to concatenate\n', ndays);

for i=1:ndays
    gfiles = dir(fullfile(days(i).name,'*.rtp'));
    [gpath, gfile, ext] = fileparts(gfiles(1).name);
    C = strsplit(gfile, '_');
    outfile = sprintf('cris2_ecmwf_csarta_clear_%s.rtp', ...
                      C{5});

   
    [h,ha,p,pa] = cat_rtp_dir(fullfile('.', days(i).name), ['*' ...
                        '.rtp']);
    nobs = length(p.rtime)

    if nobs > maxobs
        pind = randperm(nobs, maxobs);
        p2 = rtp_sub_prof(p, pind);
        p = p2;
        clear p2
    end

    rtpwrite(outfile, h,ha,p,pa);
end