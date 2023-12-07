drop table public.virus_metal_1 cascade;
create table public.virus_metal_1 as (
    select atomname,c.pdbid,chainid,resseq,residueid_ion,resname_ion,atomname_ion,atomid_ion,
           chainids_lig,resseqs_lig,residueids_lig,resnames_lig,atomnames_lig,atomids_lig,atomtypes_lig,
           resolution,header,exp_method,which_metal,coordnum_inner, qv, qc, qe, a.pdbfileid, a.bindingsiteid
    from (select pdbfileid,bindingsiteid,quality_valence as qv,quality_complete as qc,quality_experiment as qe,which_metal from neighborhood.ion_bindingsite_profiles) a
    left join (select pdbfileid,bindingsiteid,residueid_ion,resname_ion,initcap(trim(substr(atomname_ion,1,2),'-')) as atomname, trim(atomname_ion,'-') as atomname_ion,atomid_ion,protons_ion,coordnum_inner from neighborhood.ion_bindingsites) b on a.pdbfileid=b.pdbfileid and a.bindingsiteid=b.bindingsiteid
    left join (select pdbfileid,pdbid,resolution,header,exp_method from neighborhood.pdbfiles) c on a.pdbfileid=c.pdbfileid
    left join (select pdbfileid,residueid,chainid,resseq from neighborhood.residues) d on a.pdbfileid=d.pdbfileid and b.residueid_ion=d.residueid
    left join (select pdbfileid,bindingsiteid,string_to_array(conc_comma(chainid::text),', ') as chainids_lig,
                                              string_to_array(conc_comma(resseq::text),', ') as resseqs_lig,
                                              string_to_array(conc_comma(trim(resname,'-')::text),', ') as resnames_lig,
                                              string_to_array(conc_comma(residueid::text),', ') as residueids_lig,
                                              string_to_array(conc_comma(trim(atomname,'-')::text),', ') as atomnames_lig,
                                              string_to_array(conc_comma(atomid::text),', ') as atomids_lig,
                                              string_to_array(conc_comma(atomtype::text),', ') as atomtypes_lig
               from neighborhood.ion_bindingsite_ligatoms group by pdbfileid,bindingsiteid) e
               on a.pdbfileid=e.pdbfileid and a.bindingsiteid=e.bindingsiteid
);

drop table public.virus_metal cascade;
create table public.virus_metal as select * from  (
    select distinct on (atomname,pdbid,chainid,resseq,residueid_ion,resname_ion,atomname_ion,atomid_ion,
           chainids_lig,resseqs_lig,residueids_lig,resnames_lig,atomnames_lig,atomids_lig,atomtypes_lig,
           resolution,header,exp_method,which_metal,coordnum_inner, qv, qc, qe, pdbfileid, bindingsiteid)
           atomname,pdbid,chainid,resseq,residueid_ion,resname_ion,atomname_ion,atomid_ion,
           chainids_lig,resseqs_lig,residueids_lig,resnames_lig,atomnames_lig,atomids_lig,atomtypes_lig,
           resolution,header,exp_method,which_metal,coordnum_inner, qv, qc, qe, pdbfileid, bindingsiteid, chainid_virus
    from (
        select *, count(*) as chainid_cnt
        from (select *, unnest(chainids_lig) chainid_virus from public.virus_metal_1) a
        group by atomname,pdbid,chainid,resseq,residueid_ion,resname_ion,atomname_ion,atomid_ion,
           chainids_lig,resseqs_lig,residueids_lig,resnames_lig,atomnames_lig,atomids_lig,atomtypes_lig,
           resolution,header,exp_method,which_metal,coordnum_inner, qv, qc, qe, pdbfileid, bindingsiteid, chainid_virus
    ) b
    order by which_metal,coordnum_inner, qv, qc, qe desc
) a where which_metal in (3,4,11,12,13) or (which_metal>=19 and which_metal<=31) or (which_metal>=37 and which_metal<=50) or (which_metal>=55 and which_metal<=84) or which_metal>=87;

alter table virus_metal alter column exp_method type text;
update virus_metal set exp_method='XRAY' where exp_method='1';
update virus_metal set exp_method='NMR' where exp_method='2';
update virus_metal set exp_method='ELEMICRO' where exp_method='3';
update virus_metal set exp_method='SOLSCAT' where exp_method='5';

alter table public.virus_metal add column bench boolean default false;
update public.virus_metal set bench=true where qv>=0.25 and qc>=0.355 and qe>=0.5;
update public.virus_metal set bench=true where qv>=0.25 and qc>=0.355 and exp_method != 'XRAY';

alter table public.virus_metal add column bench2 boolean default false;
update public.virus_metal set bench2=true where qv>=0.5 and qc>=0.5 and qe>=0.5;
update public.virus_metal set bench2=true where qv>=0.5 and qc>=0.5 and exp_method != 'XRAY';

drop table virus_zinc;
create table virus_zinc as
        (select pdbfileid,pdbid, chainid, resseq,residueid_ion, resname_ion, atomname_ion,atomid_ion, unnest(chainids_lig) as chainid_lig, unnest(resseqs_lig) as resseq_lig,
                unnest(residueids_lig) as residueid_lig ,unnest(resnames_lig) as resname_lig, unnest(atomnames_lig) as atomname_lig,unnest(atomids_lig) as atomid_lig,unnest(atomtypes_lig) as atomtype_lig,bench,bench2 from public.virus_metal where which_metal = 30);

alter table virus_zinc alter column resseq_lig type int using resseq_lig::integer;
alter table virus_zinc alter column atomid_lig type int using atomid_lig::integer;
alter table virus_zinc alter column atomid_ion type int using atomid_ion::integer;
alter table virus_zinc alter column residueid_lig type int using residueid_lig::integer;

drop table pdb_zinc_coordinate;
create table pdb_zinc_coordinate as (
        select t1.*, t2.x as x_lig, t2.y as y_lig, t2.z as z_lig, t3.x as x_ion, t3.y as y_ion, t3.z as z_ion
        from (select a.*,
        c.residuetype from virus_zinc a
        left join neighborhood.residues c
        on a.pdbfileid = c.pdbfileid and a.residueid_lig = c.residueid) t1
        left join neighborhood.atoms t2
        on t1.pdbfileid = t2.pdbfileid and t1.residueid_lig = t2.residueid and t1.atomid_lig=t2.atomid
         left join neighborhood.atoms t3
        on t1.pdbfileid = t3.pdbfileid and t1.residueid_ion = t3.residueid and t1.atomid_ion=t3.atomid);

alter table pdb_zinc_coordinate add column site_count integer;
update pdb_zinc_coordinate t1
        set site_count = t2.site_count from (select pdbid,residueid_ion,count(distinct residueid_lig) as site_count from pdb_zinc_coordinate
                where resname_lig in ('CYS','HIS') group by pdbid,residueid_ion) t2
        where t1.pdbid=t2.pdbid and t1.residueid_ion=t2.residueid_ion;


update pdb_zinc_coordinate set bench=bench2 where site_count > 2;



----- pre site
drop table CH_site;
create table CH_site as select a.*,b.chainid as chainid_a,b.resseq as resseq_a, c.chainid as chainid_b, c.resseq as resseq_b from (
   select pdbfileid,residueid_a,residueid_b,resname_a,resname_b,atomid_a,atomid_b,atomname_a,atomname_b,atomneighbortype,neighbortype,distance as dist_ab
          from neighborhood.atomneighbors where atomneighbortype in (1111,1112,1113,1212,1213,1313) and neighbortype in (1515,1533,3315,3333) and distance>2.4
   UNION
  select pdbfileid,
          residueid_b as residueid_a, residueid_a as residueid_b,
          resname_b as resname_a,     resname_a as resname_b,
          atomid_b as atomid_a,       atomid_a as atomid_b,
          atomname_b as atomname_a,   atomname_a as atomname_b,
          atomneighbortype,neighbortype,distance as dist_ab
          from neighborhood.atomneighbors where atomneighbortype in (1111,1112,1113,1212,1213,1313) and neighbortype in (1515,1533,3315,3333) and distance>2.4
) a left join neighborhood.residues b on a.pdbfileid = b.pdbfileid and a.residueid_a=b.residueid
    left join neighborhood.residues c on a.pdbfileid = c.pdbfileid and a.residueid_b=c.residueid;

--11: side chain, sulphur atom from CYS
--13: side chain, nitrogen atom from aromatic ring (HIS)
--12: side chain, carbon atom from aromatic ring (HIS, PHE, TYR, TRP)

--15:residue_cys
--33:resiude_his


alter table CH_site add column CA_x_a real;
alter table CH_site add column CA_y_a real;
alter table CH_site add column CA_z_a real;

update CH_site t1
        set CA_x_a = t2.x from (select * from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid;
update CH_site t1
        set CA_y_a = t2.y from (select * from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid;

update CH_site t1
        set CA_z_a = t2.z from (select * from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid;

alter table CH_site add column CB_x_a real;
alter table CH_site add column CB_y_a real;
alter table CH_site add column CB_z_a real;
update CH_site t1
        set CB_x_a = t2.x from (select * from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid;
update CH_site t1
        set CB_y_a = t2.y from (select * from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid;

update CH_site t1
        set CB_z_a = t2.z from (select * from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid;

alter table CH_site add column CA_x_b real;
alter table CH_site add column CA_y_b real;
alter table CH_site add column CA_z_b real;

update CH_site t1
        set CA_x_b = t2.x from (select * from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_b = t2.residueid;
update CH_site t1
        set CA_y_b = t2.y from (select * from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_b = t2.residueid;

update CH_site t1
        set CA_z_b = t2.z from (select * from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_b = t2.residueid;

alter table CH_site add column CB_x_b real;
alter table CH_site add column CB_y_b real;
alter table CH_site add column CB_z_b real;
update CH_site t1
        set CB_x_b = t2.x from (select * from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_b = t2.residueid;
update CH_site t1
        set CB_y_b = t2.y from (select * from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_b = t2.residueid;

update CH_site t1
        set CB_z_b = t2.z from (select * from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid = t2.pdbfileid and t1.residueid_b = t2.residueid;


drop table CH_site_result_2;
create table CH_site_result_2 as select * from CH_site where residueid_a < residueid_b;


drop table CH_site_result_3;
create table CH_site_result_3 as
select * from (select z.*, y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.dist_ab as dist_bc, y.atomid_b as atomid_c, y.atomname_b as atomname_c,
                y.ca_x_b as ca_x_c, y.ca_y_b as ca_y_c, y.ca_z_b as ca_z_c,
                y.cb_x_b as cb_x_c, y.cb_y_b as cb_y_c, y.cb_z_b as cb_z_c,
                x.dist_ab as dist_ac
from CH_site z
join CH_site y on z.pdbfileid=y.pdbfileid and z.atomid_b=y.atomid_a
join CH_site x on y.pdbfileid=x.pdbfileid and y.atomid_b=x.atomid_a and x.atomid_b=z.atomid_a) t
where residueid_a < residueid_b and residueid_b < residueid_c;


drop table CH_site_result_4;
create table CH_site_result_4 as
select * from (select z.*,  y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.atomid_b as atomid_c, y.atomname_b as atomname_c,
                y.ca_x_b as ca_x_c, y.ca_y_b as ca_y_c, y.ca_z_b as ca_z_c,
                y.cb_x_b as cb_x_c, y.cb_y_b as cb_y_c, y.cb_z_b as cb_z_c,
                y.dist_ab as dist_bc,
                x.dist_ab as dist_ac, w.dist_ab as dist_ad,
                v.dist_ab as dist_bd, u.dist_ab as dist_cd,
                w.residueid_b as residueid_d,
                w.chainid_b as chainid_d, w.resseq_b as resseq_d, w.resname_b as resname_d,
                w.atomid_b as atomid_d, w.atomname_b as atomname_d,
                w.ca_x_b as ca_x_d, w.ca_y_b as ca_y_d, w.ca_z_b as ca_z_d,
                w.cb_x_b as cb_x_d, w.cb_y_b as cb_y_d, w.cb_z_b as cb_z_d
from CH_site z
join CH_site y on z.pdbfileid=y.pdbfileid and z.atomid_b=y.atomid_a
join CH_site x on y.pdbfileid=x.pdbfileid and y.atomid_b=x.atomid_a and x.atomid_b=z.atomid_a
join CH_site w on z.pdbfileid=w.pdbfileid and z.atomid_a=w.atomid_a
join CH_site v on z.pdbfileid=v.pdbfileid and z.atomid_b=v.atomid_a and w.atomid_b=v.atomid_b
join CH_site u on v.pdbfileid=u.pdbfileid and v.atomid_b=u.atomid_a and u.atomid_b=x.atomid_a) t
where residueid_a < residueid_b and residueid_b < residueid_c and residueid_c < residueid_d;

drop table ch_site_2;
create table ch_site_2 as (
        select a2.* from (select distinct a2.* from ch_site_result_2 a2
        left join ch_site_result_3 a3 on a2.pdbfileid=a3.pdbfileid and (a2.residueid_a=a3.residueid_a or a2.residueid_a=a3.residueid_b)
        where a3.pdbfileid is null) a2
        left join ch_site_result_4 a4 on a2.pdbfileid=a4.pdbfileid and (a2.residueid_a=a4.residueid_a or a2.residueid_a=a4.residueid_b or a2.residueid_a=a4.residueid_c)
        where a4.pdbfileid is null
);


drop table ch_site_3;
create table ch_site_3 as (
        select a3.* from ch_site_result_3 a3
        left join ch_site_result_4 a4 on a3.pdbfileid=a4.pdbfileid and (a3.residueid_a=a4.residueid_a or a3.residueid_a=a4.residueid_b)
        where a4.pdbfileid is NULL
);

drop table ch_site234;
create table ch_site234 as (
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
                 chainid_d,resname_d,resseq_d,residueid_d,
                 atomid_a,atomid_b,atomid_c,atomid_d,
                atomname_a,atomname_b,atomname_c,atomname_d,
               ca_x_a,ca_x_b,ca_x_c,ca_x_d,ca_y_a,ca_y_b,ca_y_c,ca_y_d,ca_z_a,ca_z_b,ca_z_c,ca_z_d,
                cb_x_a,cb_x_b,cb_x_c,cb_x_d,cb_y_a,cb_y_b,cb_y_c,cb_y_d,cb_z_a,cb_z_b,cb_z_c,cb_z_d,
                dist_ab,dist_ac,dist_bc,dist_ad,dist_bd,dist_cd
                from CH_site_result_4
    union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
                'x' as chainid_d,'x' as resname_d, -9999 as resseq_d,-9999 as residueid_d,
                atomid_a,atomid_b,atomid_c,-9999 as atomid_d,
                atomname_a,atomname_b,atomname_c,'x' as atomname_d,
               ca_x_a,ca_x_b,ca_x_c,-9999 as ca_x_d,ca_y_a,ca_y_b,ca_y_c, -9999 as ca_y_d,ca_z_a,ca_z_b,ca_z_c, -9999 as ca_z_d,
                cb_x_a,cb_x_b,cb_x_c, -9999 as cb_x_d,cb_y_a,cb_y_b,cb_y_c, -9999 as cb_y_d,cb_z_a,cb_z_b,cb_z_c, -9999 as cb_z_d,
                dist_ab,dist_ac,dist_bc,-9999 as dist_ad, -9999 as dist_bd, -9999 as dist_cd from ch_site_3
    union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                chainid_b,resname_b,resseq_b,residueid_b,
                'x' as chainid_c,'x' as resname_c, -9999 as resseq_c,-9999 as residueid_c,
                'x' as chainid_d,'x' as resname_d, -9999 as resseq_d,-9999 as residueid_d,
                atomid_a,atomid_b,-9999 as atomid_c,-9999 as atomid_d,
                atomname_a,atomname_b,'x' as atomname_c,'x' as atomname_d,
               ca_x_a,ca_x_b,-9999 as ca_x_c,-9999 as ca_x_d,ca_y_a,ca_y_b,-9999 as ca_y_c, -9999 as ca_y_d,ca_z_a,ca_z_b,-9999 as ca_z_c, -9999 as ca_z_d,
                cb_x_a,cb_x_b,-9999 as cb_x_c, -9999 as cb_x_d,cb_y_a,cb_y_b,-9999 as cb_y_c, -9999 as cb_y_d,cb_z_a,cb_z_b,-9999 as cb_z_c, -9999 as cb_z_d,
                dist_ab, -9999 as dist_ac, -9999 as dist_bc,-9999 as dist_ad, -9999 as dist_bd, -9999 as dist_cd from ch_site_2
);

drop table zinc_predict_site234;
create table zinc_predict_site234 as (
        select t1.*, t2.x as x_a, t2.y as y_a, t2.z as z_a,
        t3.x as x_b, t3.y as y_b, t3.z as z_b,
        t4.x as x_c, t4.y as y_c, t4.z as z_c,
        t5.x as x_d, t5.y as y_d, t5.z as z_d
        from (select distinct b.pdbid,a.pdbfileid,
        chainid_a, a.resseq_a, residueid_a, resname_a,a.atomid_a,a.atomname_a,
        chainid_b, a.resseq_b, residueid_b, resname_b,a.atomid_b,a.atomname_b,
        chainid_c, a.resseq_c, residueid_c, resname_c,a.atomid_c,a.atomname_c,
        chainid_d, a.resseq_d, residueid_d, resname_d,a.atomid_d,a.atomname_d,
        ca_x_a,ca_x_b,ca_x_c,ca_x_d,ca_y_a,ca_y_b,ca_y_c,ca_y_d,ca_z_a,ca_z_b,ca_z_c,ca_z_d,
        cb_x_a,cb_x_b,cb_x_c,cb_x_d,cb_y_a,cb_y_b,cb_y_c,cb_y_d,cb_z_a,cb_z_b,cb_z_c,cb_z_d,
        dist_ab
        from (select * from ch_site234)a
        left join neighborhood.pdbfiles b
        on a.pdbfileid=b.pdbfileid) t1
        left join neighborhood.atoms t2
        on t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid and t1.atomid_a = t2.atomid
         left join neighborhood.atoms t3
        on t1.pdbfileid = t3.pdbfileid and t1.residueid_b = t3.residueid and t1.atomid_b = t3.atomid
         left join neighborhood.atoms t4
        on t1.pdbfileid = t4.pdbfileid and t1.residueid_c = t4.residueid and t1.atomid_c = t4.atomid
         left join neighborhood.atoms t5
        on t1.pdbfileid = t5.pdbfileid and t1.residueid_d = t5.residueid and t1.atomid_d = t5.atomid );

alter table zinc_predict_site234 add column site_count smallint;
update zinc_predict_site234 set site_count=4 where residueid_d != -9999;
update zinc_predict_site234 set site_count=2 where residueid_c=-9999;
update zinc_predict_site234 set site_count=3 where residueid_c != -9999 and residueid_d=-9999;

delete from zinc_predict_site234 where chainid_a != chainid_b and site_count=2; 
delete from zinc_predict_site234 where atomname_a in ('-CG-') or atomname_b in ('-CG-')
      or atomname_c in ('-CG-') or atomname_d in ('-CG-') ;


alter table zinc_predict_site234 add column id integer;
update zinc_predict_site234 a set
        id = b.num from (SELECT ROW_NUMBER() over(ORDER bY pdbid,residueid_a,residueid_b,residueid_c,residueid_d DESC ) As num, 
		pdbid,residueid_a,residueid_b,residueid_c,residueid_d from zinc_predict_site234) b
        where a.pdbid=b.pdbid and a.residueid_a=b.residueid_a and a.residueid_b=b.residueid_b and a.residueid_c=b.residueid_c and a.residueid_d=b.residueid_d;

alter table zinc_predict_site234 add column resi_type text;

update zinc_predict_site234 set resi_type = 'C_C' where (resname_a = 'CYS' and resname_b = 'CYS' and resname_c = 'CYS' and (resname_d = 'CYS' or resname_d='x') );
update zinc_predict_site234 set resi_type = 'C_C' where (resname_a = 'CYS' and resname_b = 'CYS' and resname_c = 'x' and resname_d='x');

update zinc_predict_site234 set resi_type = 'H_H' where (resname_a = 'HIS' and resname_b = 'HIS' and resname_c = 'HIS' and (resname_d = 'HIS' or resname_d='x') );
update zinc_predict_site234 set resi_type = 'H_H' where (resname_a = 'HIS' and resname_b = 'HIS' and resname_c = 'x' and  resname_d='x');

update zinc_predict_site234 set resi_type = 'C_H' where resi_type is null;



drop table pre_zinc_coord_site234;
create table pre_zinc_coord_site234 as
        select distinct id,pdbid, pdbfileid, 
	resname_a, residueid_a,
	resname_b, residueid_b, 
	resname_c, residueid_c,
	resname_d, residueid_d,
        ca_x_a,ca_y_a,ca_z_a,cb_x_a,cb_y_a,cb_z_a,ca_x_b,ca_y_b,ca_z_b,cb_x_b,cb_y_b,cb_z_b,
	ca_x_c,ca_y_c,ca_z_c,cb_x_c,cb_y_c,cb_z_c,ca_x_d,ca_y_d,ca_z_d,cb_x_d,cb_y_d,cb_z_d,
        resi_type,site_count from zinc_predict_site234; 

----- lig coord
alter table pre_zinc_coord_site234 add column CD2_x float;
alter table pre_zinc_coord_site234 add column CD2_y float;
alter table pre_zinc_coord_site234 add column CD2_z float;
update pre_zinc_coord_site234 t1 set CD2_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CD2';

alter table pre_zinc_coord_site234 add column CE1_x float;
alter table pre_zinc_coord_site234 add column CE1_y float;
alter table pre_zinc_coord_site234 add column CE1_z float;
update pre_zinc_coord_site234 t1 set CE1_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CE1';

alter table pre_zinc_coord_site234 add column CG_x float;
alter table pre_zinc_coord_site234 add column CG_y float;
alter table pre_zinc_coord_site234 add column CG_z float;
update pre_zinc_coord_site234 t1 set CG_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-CG-';

alter table pre_zinc_coord_site234 add column N_x float;
alter table pre_zinc_coord_site234 add column N_y float;
alter table pre_zinc_coord_site234 add column N_z float;
update pre_zinc_coord_site234 t1 set N_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-N--';

alter table pre_zinc_coord_site234 add column ND1_x float;
alter table pre_zinc_coord_site234 add column ND1_y float;
alter table pre_zinc_coord_site234 add column ND1_z float;
update pre_zinc_coord_site234 t1 set ND1_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-ND1';

alter table pre_zinc_coord_site234 add column SG_x float;
alter table pre_zinc_coord_site234 add column SG_y float;
alter table pre_zinc_coord_site234 add column SG_z float;
update pre_zinc_coord_site234 t1 set SG_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-SG-';

alter table pre_zinc_coord_site234 add column O_x float;
alter table pre_zinc_coord_site234 add column O_y float;
alter table pre_zinc_coord_site234 add column O_z float;
update pre_zinc_coord_site234 t1 set O_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-O--';

alter table pre_zinc_coord_site234 add column NE2_x float;
alter table pre_zinc_coord_site234 add column NE2_y float;
alter table pre_zinc_coord_site234 add column NE2_z float;
update pre_zinc_coord_site234 t1 set NE2_x = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_y = t2.y from neighborhood.atoms t2
	where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-NE2';


alter table pre_zinc_coord_site234 add column C_x float;
alter table pre_zinc_coord_site234 add column C_y float;
alter table pre_zinc_coord_site234 add column C_z float;
update pre_zinc_coord_site234 t1 set C_x = t2.x from neighborhood.atoms t2
	where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_y = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_z = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t2.atomname='-C--';

------ lig1 coord
alter table pre_zinc_coord_site234 add column CD2_x1 float;
alter table pre_zinc_coord_site234 add column CD2_y1 float;
alter table pre_zinc_coord_site234 add column CD2_z1 float;
update pre_zinc_coord_site234 t1 set CD2_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CD2';

alter table pre_zinc_coord_site234 add column CE1_x1 float;
alter table pre_zinc_coord_site234 add column CE1_y1 float;
alter table pre_zinc_coord_site234 add column CE1_z1 float;
update pre_zinc_coord_site234 t1 set CE1_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CE1';

alter table pre_zinc_coord_site234 add column CG_x1 float;
alter table pre_zinc_coord_site234 add column CG_y1 float;
alter table pre_zinc_coord_site234 add column CG_z1 float;
update pre_zinc_coord_site234 t1 set CG_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-CG-';

alter table pre_zinc_coord_site234 add column N_x1 float;
alter table pre_zinc_coord_site234 add column N_y1 float;
alter table pre_zinc_coord_site234 add column N_z1 float;
update pre_zinc_coord_site234 t1 set N_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-N--';

alter table pre_zinc_coord_site234 add column ND1_x1 float;
alter table pre_zinc_coord_site234 add column ND1_y1 float;
alter table pre_zinc_coord_site234 add column ND1_z1 float;
update pre_zinc_coord_site234 t1 set ND1_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-ND1';

alter table pre_zinc_coord_site234 add column SG_x1 float;
alter table pre_zinc_coord_site234 add column SG_y1 float;
alter table pre_zinc_coord_site234 add column SG_z1 float;
update pre_zinc_coord_site234 t1 set SG_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-SG-';

alter table pre_zinc_coord_site234 add column O_x1 float;
alter table pre_zinc_coord_site234 add column O_y1 float;
alter table pre_zinc_coord_site234 add column O_z1 float;
update pre_zinc_coord_site234 t1 set O_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-O--';

alter table pre_zinc_coord_site234 add column NE2_x1 float;
alter table pre_zinc_coord_site234 add column NE2_y1 float;
alter table pre_zinc_coord_site234 add column NE2_z1 float;
update pre_zinc_coord_site234 t1 set NE2_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-NE2';


alter table pre_zinc_coord_site234 add column C_x1 float;
alter table pre_zinc_coord_site234 add column C_y1 float;
alter table pre_zinc_coord_site234 add column C_z1 float;
update pre_zinc_coord_site234 t1 set C_x1 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_y1 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_z1 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t2.atomname='-C--';

--- lig2 coord
alter table pre_zinc_coord_site234 add column CD2_x2 float;
alter table pre_zinc_coord_site234 add column CD2_y2 float;
alter table pre_zinc_coord_site234 add column CD2_z2 float;
update pre_zinc_coord_site234 t1 set CD2_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CD2';

alter table pre_zinc_coord_site234 add column CE1_x2 float;
alter table pre_zinc_coord_site234 add column CE1_y2 float;
alter table pre_zinc_coord_site234 add column CE1_z2 float;
update pre_zinc_coord_site234 t1 set CE1_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CE1';

alter table pre_zinc_coord_site234 add column CG_x2 float;
alter table pre_zinc_coord_site234 add column CG_y2 float;
alter table pre_zinc_coord_site234 add column CG_z2 float;
update pre_zinc_coord_site234 t1 set CG_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-CG-';

alter table pre_zinc_coord_site234 add column N_x2 float;
alter table pre_zinc_coord_site234 add column N_y2 float;
alter table pre_zinc_coord_site234 add column N_z2 float;
update pre_zinc_coord_site234 t1 set N_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-N--';

alter table pre_zinc_coord_site234 add column ND1_x2 float;
alter table pre_zinc_coord_site234 add column ND1_y2 float;
alter table pre_zinc_coord_site234 add column ND1_z2 float;
update pre_zinc_coord_site234 t1 set ND1_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-ND1';

alter table pre_zinc_coord_site234 add column SG_x2 float;
alter table pre_zinc_coord_site234 add column SG_y2 float;
alter table pre_zinc_coord_site234 add column SG_z2 float;
update pre_zinc_coord_site234 t1 set SG_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-SG-';

alter table pre_zinc_coord_site234 add column O_x2 float;
alter table pre_zinc_coord_site234 add column O_y2 float;
alter table pre_zinc_coord_site234 add column O_z2 float;
update pre_zinc_coord_site234 t1 set O_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-O--';

alter table pre_zinc_coord_site234 add column NE2_x2 float;
alter table pre_zinc_coord_site234 add column NE2_y2 float;
alter table pre_zinc_coord_site234 add column NE2_z2 float;
update pre_zinc_coord_site234 t1 set NE2_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-NE2';


alter table pre_zinc_coord_site234 add column C_x2 float;
alter table pre_zinc_coord_site234 add column C_y2 float;
alter table pre_zinc_coord_site234 add column C_z2 float;
update pre_zinc_coord_site234 t1 set C_x2 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_y2 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_z2 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t2.atomname='-C--';

---- lig3 coord

alter table pre_zinc_coord_site234 add column CD2_x3 float;
alter table pre_zinc_coord_site234 add column CD2_y3 float;
alter table pre_zinc_coord_site234 add column CD2_z3 float;
update pre_zinc_coord_site234 t1 set CD2_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CD2';
update pre_zinc_coord_site234 t1 set CD2_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CD2';

alter table pre_zinc_coord_site234 add column CE1_x3 float;
alter table pre_zinc_coord_site234 add column CE1_y3 float;
alter table pre_zinc_coord_site234 add column CE1_z3 float;
update pre_zinc_coord_site234 t1 set CE1_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CE1';
update pre_zinc_coord_site234 t1 set CE1_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CE1';

alter table pre_zinc_coord_site234 add column CG_x3 float;
alter table pre_zinc_coord_site234 add column CG_y3 float;
alter table pre_zinc_coord_site234 add column CG_z3 float;
update pre_zinc_coord_site234 t1 set CG_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CG-';
update pre_zinc_coord_site234 t1 set CG_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-CG-';

alter table pre_zinc_coord_site234 add column N_x3 float;
alter table pre_zinc_coord_site234 add column N_y3 float;
alter table pre_zinc_coord_site234 add column N_z3 float;
update pre_zinc_coord_site234 t1 set N_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-N--';
update pre_zinc_coord_site234 t1 set N_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-N--';

alter table pre_zinc_coord_site234 add column ND1_x3 float;
alter table pre_zinc_coord_site234 add column ND1_y3 float;
alter table pre_zinc_coord_site234 add column ND1_z3 float;
update pre_zinc_coord_site234 t1 set ND1_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-ND1';
update pre_zinc_coord_site234 t1 set ND1_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-ND1';

alter table pre_zinc_coord_site234 add column SG_x3 float;
alter table pre_zinc_coord_site234 add column SG_y3 float;
alter table pre_zinc_coord_site234 add column SG_z3 float;
update pre_zinc_coord_site234 t1 set SG_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-SG-';
update pre_zinc_coord_site234 t1 set SG_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-SG-';

alter table pre_zinc_coord_site234 add column O_x3 float;
alter table pre_zinc_coord_site234 add column O_y3 float;
alter table pre_zinc_coord_site234 add column O_z3 float;
update pre_zinc_coord_site234 t1 set O_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-O--';
update pre_zinc_coord_site234 t1 set O_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-O--';

alter table pre_zinc_coord_site234 add column NE2_x3 float;
alter table pre_zinc_coord_site234 add column NE2_y3 float;
alter table pre_zinc_coord_site234 add column NE2_z3 float;
update pre_zinc_coord_site234 t1 set NE2_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-NE2';
update pre_zinc_coord_site234 t1 set NE2_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-NE2';


alter table pre_zinc_coord_site234 add column C_x3 float;
alter table pre_zinc_coord_site234 add column C_y3 float;
alter table pre_zinc_coord_site234 add column C_z3 float;
update pre_zinc_coord_site234 t1 set C_x3 = t2.x from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_y3 = t2.y from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-C--';

update pre_zinc_coord_site234 t1 set C_z3 = t2.z from neighborhood.atoms t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t2.atomname='-C--';


alter table pre_zinc_coord_site234 add column ab_type text;
update pre_zinc_coord_site234 set ab_type='H_H' where resname_a='HIS' and resname_b='HIS';
update pre_zinc_coord_site234 set ab_type='C_C' where resname_a='CYS' and resname_b='CYS';
update pre_zinc_coord_site234 set ab_type='C_H' where (resname_a='CYS' and resname_b='HIS') or (resname_a='HIS' and resname_b='CYS');

alter table pre_zinc_coord_site234 add column ac_type text;
update pre_zinc_coord_site234 set ac_type='H_H' where resname_a='HIS' and resname_c='HIS';
update pre_zinc_coord_site234 set ac_type='C_C' where resname_a='CYS' and resname_c='CYS';
update pre_zinc_coord_site234 set ac_type='C_H' where (resname_a='CYS' and resname_c='HIS') or (resname_a='HIS' and resname_c='CYS');

alter table pre_zinc_coord_site234 add column bc_type text;
update pre_zinc_coord_site234 set bc_type='H_H' where resname_b='HIS' and resname_c='HIS';
update pre_zinc_coord_site234 set bc_type='C_C' where resname_b='CYS' and resname_c='CYS';
update pre_zinc_coord_site234 set bc_type='C_H' where (resname_b='CYS' and resname_c='HIS') or (resname_b='HIS' and resname_c='CYS');

