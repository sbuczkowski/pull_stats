rtp_addpaths
year = 2017;
basedir = '/asl/rtp_lustre/rtp_airibrad_v5/2017/random';

rtpfiles = dir(fullfile(basedir, 'era_*.rtp'));
fprintf(1, '> Found %d files to check and normalize\n', ...
        length(rtpfiles));

for ifile = 1:length(rtpfiles)
    [h,ha,p,pa] = rtpread(fullfile(rtpfiles(ifile).folder, ...
                          rtpfiles(ifile).name));

    fprintf(1, '> Reading and verifying %s\n', rtpfiles(ifile).name);
    % check for existence of old fieldnames. If they exist, copy
    % field to new fieldname (assumes that not finding old
    % fieldnames implies existence of new fieldnames). Once copied,
    % remove the old field
    if (isfield(p, 'rcalc') | isfield(p, 'sarta_rclearcalc'))
        if isfield(p, 'rcalc')
            fprintf(1, '>> Found rcalc -> rcld\n');
            p.rcld = p.rcalc;
            p = rmfield(p, 'rcalc');
        end
        if isfield(p, 'sarta_rclearcalc')
            fprintf(1, '>> Found sarta_rclearcalc -> rclr\n');
            p.rclr = p.sarta_rclearcalc;
            p = rmfield(p, 'sarta_rclearcalc');
        end
        fprintf(1, '>> Structure conforms to new standard\n');
        
        rtpwrite(fullfile(rtpfiles(ifile).folder, sprintf('%s.1', ...
                                                          rtpfiles(ifile).name)), h,ha,p,pa);
    end
    fprintf(1, '>> File conforms to new standard\n');

end

fprintf(1, '> Done\n');
