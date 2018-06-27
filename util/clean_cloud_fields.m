function clean_cloud_fields(p)

        % Remove 'clouds' that have only partial defining
        % characteristics (effects 0.2% of obs in test rtp data)
        % cfrac=cngwat=0 but cpsize/cprtop/cprbot ~= 0
        a = find(p.cfrac==0 & p.cngwat==0 & p.cprtop>0);
        p.ctype(a)=-1;
        a = find(p.cfrac2==0 & p.cngwat2==0 & p.cprtop2>0);
        p.ctype2(a)=-1;

        % replace 'bad data' placeholders in cprtop/cprbot (-9999
        % -> NaN, esesntially except for 'clouds' from above)
        a = find(p.ctype==-1);
        p.cprtop(a)=NaN; p.cprbot(a)=NaN;
        a = find(p.ctype2==-1);
        p.cprtop2(a)=NaN; p.cprbot2(a)=NaN;

        % massage situations where we have two of same cloud type
        % in a single obs i.e. p.ctype=p.ctype2 (average to a
        % single cloud)
        dblwat = find(p.ctype==101 & p.ctype2==101);
        p.cfrac2(dblwat) = p.cfrac(dblwat) + p.cfrac2(dblwat);
        p.cpsize2(dblwat) = (p.cpsize(dblwat) + ...
                             p.cpsize2(dblwat))/2;
        p.cprtop2(dblwat) = (p.cprtop(dblwat) + ...
                             p.cprtop2(dblwat))/2;
        p.cprbot2(dblwat) = (p.cprbot(dblwat) + ...
                             p.cprbot2(dblwat))/2;
        p.cngwat2(dblwat) = (p.cfrac(dblwat).*p.cngwat(dblwat) + ...
                             p.cfrac2(dblwat).*p.cngwat2(dblwat))./ ...
            p.cfrac2(dblwat);
        p.cfrac(dblwat)=0; p.cpsize(dblwat)=0;
        p.cngwat(dblwat)= 0; p.cprtop(dblwat)=NaN;
        p.cprbot(dblwat)=NaN; p.ctype(dblwat)=-1;
        
        dblice = find(p.ctype==201 & p.ctype2==201);
        p.cfrac(dblice) = p.cfrac(dblice) + p.cfrac2(dblice);
        p.cpsize(dblice) = (p.cpsize(dblice) + ...
                            p.cpsize2(dblice))/2;
        p.cprtop(dblice) = (p.cprtop(dblice) + ...
                            p.cprtop2(dblice))/2;
        p.cprbot(dblice) = (p.cprbot(dblice) + ...
                            p.cprbot2(dblice))/2;
        p.cngwat(dblice) = (p.cfrac(dblice).*p.cngwat(dblice) + ...
                            p.cfrac2(dblice).*p.cngwat2(dblice))./ ...
            p.cfrac2(dblice);
        p.cfrac2(dblice)=0; p.cpsize2(dblice)=0;
        p.cngwat2(dblice)= 0; p.cprtop2(dblice)=NaN;
        p.cprbot2(dblice)=NaN; p.ctype2(dblice)=-1;
        
        % at this point cloud fields should now be set so that
        % ctype/cpsize/cngwat and similar contain only ice clouds
        % (or none) and ctype2/cpsize2/cngwat2 and similar
        % contain only water clouds (or none)
        trace.NOTE_ON_CLOUD_VARS = ['ctype_mean and similar are average ' ...
                            'ICE cloud parameters, ctype2_mean and similar ' ...
                            'are average WATER cloud parameters'];

        % pull out cloud parameters for averaging
        % build logical array for each cloud type and pull out
        % only values associated with a cloud for averaging
        lIce = (p.ctype==201 & p.cfrac>0 & p.cngwat>0 & p.cpsize>0 ...
                & p.cprbot>0 & p.cprtop>0);

        ice_ind = find(lIce);
        count_ice(iday, ilat) = length(ice_ind);
        ctype_mean(iday, ilat) = nanmean(p.ctype(ice_ind));
        cngwat_mean(iday, ilat) = nanmean(p.cngwat(ice_ind));
        cpsize_mean(iday, ilat) = nanmean(p.cpsize(ice_ind));
        cprbot_mean(iday, ilat) = nanmean(p.cprbot(ice_ind));
        cprtop_mean(iday, ilat) = nanmean(p.cprtop(ice_ind));
        
        lWat = (p.ctype2==101 & p.cfrac2>0 & p.cngwat2>0 & p.cpsize2>0 ...
                & p.cprbot2>0 & p.cprtop2>0);

        wat_ind = find(lWat);
        count_water(iday, ilat) = length(wat_ind); 
        ctype2_mean(iday, ilat) = nanmean(p.ctype2(wat_ind));
        cngwat2_mean(iday, ilat) = nanmean(p.cngwat2(wat_ind));
        cpsize2_mean(iday, ilat) = nanmean(p.cpsize2(wat_ind));
        cprbot2_mean(iday, ilat) = nanmean(p.cprbot2(wat_ind));
        cprtop2_mean(iday, ilat) = nanmean(p.cprtop2(wat_ind));

        % cloud fraction gets averaged over ALL obs
        cfrac_mean(iday, ilat) = nanmean(p.cfrac);
        cfrac2_mean(iday, ilat) = nanmean(p.cfrac2);
        cfrac12_mean(iday, ilat) = nanmean(p.cfrac12);

        