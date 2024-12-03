--10 - side chain, oxygen atom from carboxylic groups of ASP & GLU
--11 - side chain, sulphur atom from CYS
--12 - side chain, carbon atom from aromatic ring (HIS, PHE, TYR, TRP)  (CE1, CD2)
--13 - side chain, nitrogen atom from aromatic ring (HIS)               (ND1, NE2)

--define RESIDUE_CYS             15
--define RESIDUE_ASP             21
--define RESIDUE_GLU             25
--define RESIDUE_HIS             33

---- 1. The existing transition metal binding sites were obtained.
drop table all_metal_1 cascade;
create table all_metal_1 as (
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

drop table all_metal cascade;
create table all_metal as select * from  (
    select distinct on (atomname,pdbid,chainid,resseq,residueid_ion,resname_ion,atomname_ion,atomid_ion,
           chainids_lig,resseqs_lig,residueids_lig,resnames_lig,atomnames_lig,atomids_lig,atomtypes_lig,
           resolution,header,exp_method,which_metal,coordnum_inner, qv, qc, qe, pdbfileid, bindingsiteid)
           atomname,pdbid,chainid,resseq,residueid_ion,resname_ion,atomname_ion,atomid_ion,
           chainids_lig,resseqs_lig,residueids_lig,resnames_lig,atomnames_lig,atomids_lig,atomtypes_lig,
           resolution,header,exp_method,which_metal,coordnum_inner, qv, qc, qe, pdbfileid, bindingsiteid, chainid_virus
    from (
        select *, count(*) as chainid_cnt
        from (select *, unnest(chainids_lig) chainid_virus from all_metal_1) a
        group by atomname,pdbid,chainid,resseq,residueid_ion,resname_ion,atomname_ion,atomid_ion,
           chainids_lig,resseqs_lig,residueids_lig,resnames_lig,atomnames_lig,atomids_lig,atomtypes_lig,
           resolution,header,exp_method,which_metal,coordnum_inner, qv, qc, qe, pdbfileid, bindingsiteid, chainid_virus
    ) b
    order by which_metal,coordnum_inner, qv, qc, qe desc
) a where which_metal in (3,4,11,12,13) or (which_metal>=19 and which_metal<=31) or (which_metal>=37 and which_metal<=50) or (which_metal>=55 and which_metal<=84) or which_metal>=87;

alter table all_metal alter column exp_method type text;
update all_metal set exp_method='XRAY' where exp_method='1';
update all_metal set exp_method='NMR' where exp_method='2';
update all_metal set exp_method='ELEMICRO' where exp_method='3';
update all_metal set exp_method='SOLSCAT' where exp_method='5';

alter table all_metal add column bench boolean default false;
update all_metal set bench=true where qv>=0.25 and qc>=0.355 and qe>=0.5;
update all_metal set bench=true where qv>=0.25 and qc>=0.355 and exp_method != 'XRAY';

alter table all_metal add column bench2 boolean default false;
update all_metal set bench2=true where qv>=0.5 and qc>=0.5 and qe>=0.5;
update all_metal set bench2=true where qv>=0.5 and qc>=0.5 and exp_method != 'XRAY';

drop table transition_metal;
create table transition_metal as select a.*,b.occupancy_ion,b.occupancy_env_avg from
        (select pdbfileid,pdbid, chainid, resseq,residueid_ion, resname_ion, atomname_ion,atomid_ion, unnest(chainids_lig) as chainid_lig, unnest(resseqs_lig) as resseq_lig,
                unnest(residueids_lig) as residueid_lig ,unnest(resnames_lig) as resname_lig, unnest(atomnames_lig) as atomname_lig,unnest(atomids_lig) as atomid_lig,unnest(atomtypes_lig) as atomtype_lig,qv,qc,qe,bench,bench2,which_metal from all_metal
                where which_metal in (25,26,27,28,29,30)) a
	left join neighborhood.ion_bindingsites b 
	on a.pdbfileid=b.pdbfileid and a.residueid_ion=b.residueid_ion;

alter table transition_metal alter column resseq_lig type int using resseq_lig::integer;
alter table transition_metal alter column atomid_lig type int using atomid_lig::integer;
alter table transition_metal alter column atomid_ion type int using atomid_ion::integer;
alter table transition_metal alter column residueid_lig type int using residueid_lig::integer;


drop table transition_metal_coord;
create table transition_metal_coord as select t1.*,t3.x as x_ion, t3.y as y_ion, t3.z as z_ion from
         (select distinct a.pdbid,a.pdbfileid,a.chainid,a.resseq,a.residueid_ion,a.which_metal,a.atomname_ion,a.resname_ion,a.atomid_ion,a.chainid_lig,a.resseq_lig,a.residueid_lig,a.resname_lig,a.atomid_lig,a.bench,a.bench2,occupancy_ion,
        c.residuetype from (select * from transition_metal) a 
        left join neighborhood.residues c
        on a.pdbfileid = c.pdbfileid and a.residueid_lig = c.residueid) t1
        left join neighborhood.atoms t3
        on t1.pdbfileid = t3.pdbfileid and t1.residueid_ion = t3.residueid and t1.atomid_ion=t3.atomid;

alter table transition_metal_coord add column aa_count integer default 0;
update transition_metal_coord t1
        set aa_count = t2.site_count from (select pdbid,residueid_ion,count(distinct residueid_lig) as site_count from transition_metal_coord where residuetype <= 41 group by pdbid,residueid_ion) t2
        where t1.pdbid=t2.pdbid and t1.residueid_ion=t2.residueid_ion ;

update transition_metal_coord set bench=bench2 where aa_count > 2;

alter table transition_metal_coord add column id integer;
update transition_metal_coord a set
        id = b.num from (SELECT ROW_NUMBER() over(ORDER bY pdbid,residueid_ion DESC ) As num,
                pdbid,residueid_ion from transition_metal_coord) b
        where a.pdbid=b.pdbid and a.residueid_ion=b.residueid_ion;


---------2. Prediction of candidate metal binding sites (CH)

drop table CH_site;
create table CH_site as select a.*,b.chainid as chainid_a,b.resseq as resseq_a, c.chainid as chainid_b, c.resseq as resseq_b from (
   select pdbfileid,residueid_a,residueid_b,resname_a,resname_b,atomid_a,atomid_b,atomname_a,atomname_b,atomneighbortype,neighbortype,distance as dist_ab, contact_flag as contact_ab
          from neighborhood.atomneighbors where atomneighbortype in (1111,1112,1113,1212,1213,1313) and neighbortype in (1515,1533,3315,3333) and distance>2.4
   UNION
  select pdbfileid,
          residueid_b as residueid_a, residueid_a as residueid_b,
          resname_b as resname_a,     resname_a as resname_b,
          atomid_b as atomid_a,       atomid_a as atomid_b,
          atomname_b as atomname_a,   atomname_a as atomname_b,
          atomneighbortype,neighbortype,distance as dist_ab, contact_flag as contact_ab
          from neighborhood.atomneighbors where atomneighbortype in (1111,1112,1113,1212,1213,1313) and neighbortype in (1515,1533,3315,3333) and distance>2.4
) a left join neighborhood.residues b on a.pdbfileid = b.pdbfileid and a.residueid_a=b.residueid
    left join neighborhood.residues c on a.pdbfileid = c.pdbfileid and a.residueid_b=c.residueid;


drop table CH_site_result_2;
create table CH_site_result_2 as select * from CH_site where residueid_a < residueid_b;


drop table CH_site_result_3;
create table CH_site_result_3 as
select * from (select z.*, y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.atomid_b as atomid_c, y.atomname_b as atomname_c,
                y.contact_ab as contact_bc, x.contact_ab as contact_ca
from CH_site z
join CH_site y on z.pdbfileid=y.pdbfileid and z.atomid_b=y.atomid_a
join CH_site x on y.pdbfileid=x.pdbfileid and y.atomid_b=x.atomid_a and x.atomid_b=z.atomid_a) t
where residueid_a < residueid_b and residueid_b < residueid_c;


drop table CH_site_result_4;
create table CH_site_result_4 as
select * from (select z.*,  y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.atomid_b as atomid_c, y.atomname_b as atomname_c,
                w.residueid_b as residueid_d,
                w.chainid_b as chainid_d, w.resseq_b as resseq_d, w.resname_b as resname_d,
                w.atomid_b as atomid_d, w.atomname_b as atomname_d,
                y.contact_ab as contact_bc,
                x.contact_ab as contact_ca, w.contact_ab as contact_ad,
                v.contact_ab as contact_bd, u.contact_ab as contact_cd
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
create table ch_site234 as
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
                 chainid_d,resname_d,resseq_d,residueid_d,
                 atomid_a,atomid_b,atomid_c,atomid_d,
                atomname_a,atomname_b,atomname_c,atomname_d,
                contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd, 
		4 as site_count
                from CH_site_result_4
        union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
                 'NaN' as chainid_d,'NaN' as resname_d,-9999 as resseq_d,-9999 as residueid_d,
                 atomid_a,atomid_b,atomid_c,-9999 as atomid_d,
                atomname_a,atomname_b,atomname_c,'NaN' as atomname_d,
                contact_ab,contact_bc,contact_ca,cast(1111 as bit) as contact_ad,cast(1111 as bit) as contact_bd,cast(1111 as bit) as contact_cd,
		3 as site_count
                from CH_site_3
        union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 'NaN' as chainid_c,'NaN' as resname_c,-9999 as resseq_c,-9999 as residueid_c,
                 'NaN' as chainid_d,'NaN' as resname_d,-9999 as resseq_d,-9999 as residueid_d,
                 atomid_a,atomid_b,-9999 as atomid_c,-9999 as atomid_d,
                atomname_a,atomname_b,'NaN' as atomname_c,'NaN' as atomname_d,
                contact_ab,cast(1111 as bit) as contact_bc,cast(1111 as bit) as contact_ca,cast(1111 as bit) as contact_ad,cast(1111 as bit) as contact_bd,cast(1111 as bit) as contact_cd,
		2 as site_count
                from CH_site_2;



drop table zncu_predict_sites;
create table zncu_predict_sites as
        select t1.*, t2.x as x_a, t2.y as y_a, t2.z as z_a,
        t3.x as x_b, t3.y as y_b, t3.z as z_b,
        t4.x as x_c, t4.y as y_c, t4.z as z_c,
        t5.x as x_d, t5.y as y_d, t5.z as z_d
        from (select distinct b.pdbid,a.* from ch_site234 a
        left join neighborhood.pdbfiles b
        on a.pdbfileid=b.pdbfileid) t1
        left join neighborhood.atoms t2
        on t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid and t1.atomid_a = t2.atomid
         left join neighborhood.atoms t3
        on t1.pdbfileid = t3.pdbfileid and t1.residueid_b = t3.residueid and t1.atomid_b = t3.atomid
         left join neighborhood.atoms t4
        on t1.pdbfileid = t4.pdbfileid and t1.residueid_c = t4.residueid and t1.atomid_c = t4.atomid
         left join neighborhood.atoms t5
        on t1.pdbfileid = t5.pdbfileid and t1.residueid_d = t5.residueid and t1.atomid_d = t5.atomid;


----delete from zncu_predict_sites where chainid_a != chainid_b and site_count=2;
delete from zncu_predict_sites where atomname_a in ('-CG-') or atomname_b in ('-CG-')
      or atomname_c in ('-CG-') or atomname_d in ('-CG-') ;


alter table zncu_predict_sites add column id integer;
update zncu_predict_sites a set
        id = b.num from (SELECT ROW_NUMBER() over(ORDER bY pdbid,residueid_a,residueid_b,residueid_c,residueid_d DESC ) As num,
                pdbid,residueid_a,residueid_b,residueid_c,residueid_d from zncu_predict_sites) b
        where a.pdbid=b.pdbid and a.residueid_a=b.residueid_a and a.residueid_b=b.residueid_b and a.residueid_c=b.residueid_c and a.residueid_d=b.residueid_d ;

alter table zncu_predict_sites add column resi_type text;

update zncu_predict_sites set resi_type = 'C_C' where
        (resname_a in ('CYS') and resname_b in ('CYS') and site_count=2) or
        (resname_a in ('CYS') and resname_b in ('CYS') and resname_c in ('CYS') and site_count=3) or
        (resname_a in ('CYS') and resname_b in ('CYS') and resname_c in ('CYS') and resname_d in ('CYS') and site_count=4);

update zncu_predict_sites set resi_type = 'H_H' where
        (resname_a in ('HIS') and resname_b in ('HIS') and site_count=2) or
        (resname_a in ('HIS') and resname_b in ('HIS') and resname_c in ('HIS') and site_count=3) or
        (resname_a in ('HIS') and resname_b in ('HIS') and resname_c in ('HIS') and resname_d in ('HIS') and site_count=4); 

update zncu_predict_sites set resi_type = 'C_H' where resi_type is null;

------ 3. Add atomic coordinate column (CH)

drop table pre_zncu_coord_sites;
create table pre_zncu_coord_sites as
        select distinct id,pdbid, pdbfileid, 
	resname_a, residueid_a,
	resname_b, residueid_b, 
	resname_c, residueid_c,
	resname_d, residueid_d,
	contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd,
        resi_type,site_count from zncu_predict_sites; 

----- CD2 coord
alter table pre_zncu_coord_sites add column CD2_coords_a text;
update pre_zncu_coord_sites t1 set CD2_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column CD2_coords_b text;
update pre_zncu_coord_sites t1 set CD2_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column CD2_coords_c text;
update pre_zncu_coord_sites t1 set CD2_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column CD2_coords_d text;
update pre_zncu_coord_sites t1 set CD2_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------ CE1 coord
alter table pre_zncu_coord_sites add column CE1_coords_a text;
update pre_zncu_coord_sites t1 set CE1_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column CE1_coords_b text;
update pre_zncu_coord_sites t1 set CE1_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column CE1_coords_c text;
update pre_zncu_coord_sites t1 set CE1_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column CE1_coords_d text;
update pre_zncu_coord_sites t1 set CE1_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- N coord
alter table pre_zncu_coord_sites add column N_coords_a text;
update pre_zncu_coord_sites t1 set N_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column N_coords_b text;
update pre_zncu_coord_sites t1 set N_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column N_coords_c text;
update pre_zncu_coord_sites t1 set N_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column N_coords_d text;
update pre_zncu_coord_sites t1 set N_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- ND1 coord
alter table pre_zncu_coord_sites add column ND1_coords_a text;
update pre_zncu_coord_sites t1 set ND1_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column ND1_coords_b text;
update pre_zncu_coord_sites t1 set ND1_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column ND1_coords_c text;
update pre_zncu_coord_sites t1 set ND1_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column ND1_coords_d text;
update pre_zncu_coord_sites t1 set ND1_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- SG coord
alter table pre_zncu_coord_sites add column SG_coords_a text;
update pre_zncu_coord_sites t1 set SG_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column SG_coords_b text;
update pre_zncu_coord_sites t1 set SG_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column SG_coords_c text;
update pre_zncu_coord_sites t1 set SG_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column SG_coords_d text;
update pre_zncu_coord_sites t1 set SG_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- O coord
alter table pre_zncu_coord_sites add column O_coords_a text;
update pre_zncu_coord_sites t1 set O_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column O_coords_b text;
update pre_zncu_coord_sites t1 set O_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column O_coords_c text;
update pre_zncu_coord_sites t1 set O_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column O_coords_d text;
update pre_zncu_coord_sites t1 set O_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- NE2 coord
alter table pre_zncu_coord_sites add column NE2_coords_a text;
update pre_zncu_coord_sites t1 set NE2_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column NE2_coords_b text;
update pre_zncu_coord_sites t1 set NE2_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column NE2_coords_c text;
update pre_zncu_coord_sites t1 set NE2_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column NE2_coords_d text;
update pre_zncu_coord_sites t1 set NE2_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

-------- C coord
alter table pre_zncu_coord_sites add column C_coords_a text;
update pre_zncu_coord_sites t1 set C_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column C_coords_b text;
update pre_zncu_coord_sites t1 set C_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column C_coords_c text;
update pre_zncu_coord_sites t1 set C_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column C_coords_d text;
update pre_zncu_coord_sites t1 set C_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- CA coord
alter table pre_zncu_coord_sites add column CA_coords_a text;
update pre_zncu_coord_sites t1 set CA_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column CA_coords_b text;
update pre_zncu_coord_sites t1 set CA_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column CA_coords_c text;
update pre_zncu_coord_sites t1 set CA_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column CA_coords_d text;
update pre_zncu_coord_sites t1 set CA_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- CB coord
alter table pre_zncu_coord_sites add column CB_coords_a text;
update pre_zncu_coord_sites t1 set CB_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column CB_coords_b text;
update pre_zncu_coord_sites t1 set CB_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column CB_coords_c text;
update pre_zncu_coord_sites t1 set CB_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column CB_coords_d text;
update pre_zncu_coord_sites t1 set CB_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

------- CG coord
alter table pre_zncu_coord_sites add column CG_coords_a text;
update pre_zncu_coord_sites t1 set CG_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_zncu_coord_sites add column CG_coords_b text;
update pre_zncu_coord_sites t1 set CG_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_zncu_coord_sites add column CG_coords_c text;
update pre_zncu_coord_sites t1 set CG_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_zncu_coord_sites add column CG_coords_d text;
update pre_zncu_coord_sites t1 set CG_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

--------- 4. Prediction of candidate metal binding sites (EDH)

drop table EDH_site;
create table EDH_site as select a.*,b.chainid as chainid_a,b.resseq as resseq_a, c.chainid as chainid_b, c.resseq as resseq_b from (
   select pdbfileid,residueid_a,residueid_b,resname_a,resname_b,atomid_a,atomid_b,atomname_a,atomname_b,atomneighbortype,neighbortype,distance as dist_ab, contact_flag as contact_ab
          from neighborhood.atomneighbors where atomneighbortype in (1010,1012,1013,1212,1213,1313) and neighbortype in (2121,2133,3321,3333,2125,2521,2525,2533,3325) and distance>2.4
   UNION
  select pdbfileid,
          residueid_b as residueid_a, residueid_a as residueid_b,
          resname_b as resname_a,     resname_a as resname_b,
          atomid_b as atomid_a,       atomid_a as atomid_b,
          atomname_b as atomname_a,   atomname_a as atomname_b,
          atomneighbortype,neighbortype,distance as dist_ab,contact_flag as contact_ab
          from neighborhood.atomneighbors where atomneighbortype in (1010,1012,1013,1212,1213,1313) and neighbortype in (2121,2133,3321,3333,2125,2521,2525,2533,3325) and distance>2.4
) a left join neighborhood.residues b on a.pdbfileid = b.pdbfileid and a.residueid_a=b.residueid
    left join neighborhood.residues c on a.pdbfileid = c.pdbfileid and a.residueid_b=c.residueid;

drop table EDH_site_result_2;
create table EDH_site_result_2 as select * from EDH_site where residueid_a < residueid_b;

drop table EDH_site_result_3;
create table EDH_site_result_3 as
select * from (select z.*, y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.atomid_b as atomid_c, y.atomname_b as atomname_c,
		y.contact_ab as contact_bc, x.contact_ab as contact_ca
from EDH_site z
join EDH_site y on z.pdbfileid=y.pdbfileid and z.atomid_b=y.atomid_a
join EDH_site x on y.pdbfileid=x.pdbfileid and y.atomid_b=x.atomid_a and x.atomid_b=z.atomid_a) t
where residueid_a < residueid_b and residueid_b < residueid_c;

drop table EDH_site_result_4;
create table EDH_site_result_4 as
select * from (select z.*,  y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.atomid_b as atomid_c, y.atomname_b as atomname_c,
                w.residueid_b as residueid_d,
                w.chainid_b as chainid_d, w.resseq_b as resseq_d, w.resname_b as resname_d,
                w.atomid_b as atomid_d, w.atomname_b as atomname_d,
		y.contact_ab as contact_bc,
                x.contact_ab as contact_ca, w.contact_ab as contact_ad,
                v.contact_ab as contact_bd, u.contact_ab as contact_cd
from EDH_site z
join EDH_site y on z.pdbfileid=y.pdbfileid and z.atomid_b=y.atomid_a
join EDH_site x on y.pdbfileid=x.pdbfileid and y.atomid_b=x.atomid_a and x.atomid_b=z.atomid_a
join EDH_site w on z.pdbfileid=w.pdbfileid and z.atomid_a=w.atomid_a
join EDH_site v on z.pdbfileid=v.pdbfileid and z.atomid_b=v.atomid_a and w.atomid_b=v.atomid_b
join EDH_site u on v.pdbfileid=u.pdbfileid and v.atomid_b=u.atomid_a and u.atomid_b=x.atomid_a) t
where residueid_a < residueid_b and residueid_b < residueid_c and residueid_c < residueid_d;

drop table EDH_site_result_5;
create table EDH_site_result_5 as
select * from (select z.*,  y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.atomid_b as atomid_c, y.atomname_b as atomname_c,
                w.residueid_b as residueid_d,
                w.chainid_b as chainid_d, w.resseq_b as resseq_d, w.resname_b as resname_d,
                w.atomid_b as atomid_d, w.atomname_b as atomname_d,
		t.residueid_b as residueid_e,
                t.chainid_b as chainid_e, t.resseq_b as resseq_e, t.resname_b as resname_e,
                t.atomid_b as atomid_e, t.atomname_b as atomname_e,
		y.contact_ab as contact_bc,
                x.contact_ab as contact_ca, w.contact_ab as contact_ad,
                v.contact_ab as contact_bd, u.contact_ab as contact_cd,
                t.contact_ab as contact_ae, s.contact_ab as contact_be, q.contact_ab as contact_ce, q.contact_ab as contact_de
from EDH_site z
join EDH_site y on z.pdbfileid=y.pdbfileid and z.atomid_b=y.atomid_a
join EDH_site x on y.pdbfileid=x.pdbfileid and y.atomid_b=x.atomid_a and x.atomid_b=z.atomid_a
join EDH_site w on z.pdbfileid=w.pdbfileid and z.atomid_a=w.atomid_a
join EDH_site v on z.pdbfileid=v.pdbfileid and z.atomid_b=v.atomid_a and w.atomid_b=v.atomid_b
join EDH_site u on v.pdbfileid=u.pdbfileid and v.atomid_b=u.atomid_a and u.atomid_b=x.atomid_a
join EDH_site t on z.pdbfileid=t.pdbfileid and z.atomid_a=t.atomid_a
join EDH_site s on z.pdbfileid=s.pdbfileid and z.atomid_b=s.atomid_a and s.atomid_b=t.atomid_b
join EDH_site r on t.pdbfileid=r.pdbfileid and t.atomid_b=r.atomid_a and r.atomid_b=y.atomid_b
join EDH_site q on t.pdbfileid=q.pdbfileid and t.atomid_b=q.atomid_a and q.atomid_b=w.atomid_b
) t
where residueid_a < residueid_b and residueid_b < residueid_c and residueid_c < residueid_d and residueid_d < residueid_e;

drop table EDH_site_result_6;
create table EDH_site_result_6 as
select * from (select z.*,  y.residueid_b as residueid_c,
                y.chainid_b as chainid_c, y.resseq_b as resseq_c, y.resname_b as resname_c,
                y.atomid_b as atomid_c, y.atomname_b as atomname_c,
                w.residueid_b as residueid_d,
                w.chainid_b as chainid_d, w.resseq_b as resseq_d, w.resname_b as resname_d,
                w.atomid_b as atomid_d, w.atomname_b as atomname_d,
                t.residueid_b as residueid_e,
                t.chainid_b as chainid_e, t.resseq_b as resseq_e, t.resname_b as resname_e,
                t.atomid_b as atomid_e, t.atomname_b as atomname_e,
		p.residueid_b as residueid_f,
                p.chainid_b as chainid_f, p.resseq_b as resseq_f, p.resname_b as resname_f,
                p.atomid_b as atomid_f, p.atomname_b as atomname_f,
		y.contact_ab as contact_bc,
                x.contact_ab as contact_ca, w.contact_ab as contact_ad,
                v.contact_ab as contact_bd, u.contact_ab as contact_cd,
                t.contact_ab as contact_ae, s.contact_ab as contact_be, q.contact_ab as contact_ce, q.contact_ab as contact_de,
                o.contact_ab as contact_af, n.contact_ab as contact_bf, m.contact_ab as contact_cf, l.contact_ab as contact_df, p.contact_ab as contact_ef
from EDH_site z
join EDH_site y on z.pdbfileid=y.pdbfileid and z.atomid_b=y.atomid_a
join EDH_site x on y.pdbfileid=x.pdbfileid and y.atomid_b=x.atomid_a and x.atomid_b=z.atomid_a
join EDH_site w on z.pdbfileid=w.pdbfileid and z.atomid_a=w.atomid_a
join EDH_site v on z.pdbfileid=v.pdbfileid and z.atomid_b=v.atomid_a and w.atomid_b=v.atomid_b
join EDH_site u on v.pdbfileid=u.pdbfileid and v.atomid_b=u.atomid_a and u.atomid_b=x.atomid_a
join EDH_site t on z.pdbfileid=t.pdbfileid and z.atomid_a=t.atomid_a
join EDH_site s on z.pdbfileid=s.pdbfileid and z.atomid_b=s.atomid_a and s.atomid_b=t.atomid_b
join EDH_site r on t.pdbfileid=r.pdbfileid and t.atomid_b=r.atomid_a and r.atomid_b=y.atomid_b
join EDH_site q on t.pdbfileid=q.pdbfileid and t.atomid_b=q.atomid_a and q.atomid_b=w.atomid_b
join EDH_site p on t.pdbfileid=p.pdbfileid and t.atomid_b=p.atomid_a
join EDH_site o on z.pdbfileid=o.pdbfileid and z.atomid_a=o.atomid_a and o.atomid_b=p.atomid_b
join EDH_site n on z.pdbfileid=n.pdbfileid and z.atomid_b=n.atomid_a and n.atomid_b=p.atomid_b
join EDH_site m on p.pdbfileid=m.pdbfileid and p.atomid_b=m.atomid_a and m.atomid_b=y.atomid_b
join EDH_site l on p.pdbfileid=l.pdbfileid and p.atomid_b=l.atomid_a and l.atomid_b=w.atomid_b
) t
where residueid_a < residueid_b and residueid_b < residueid_c and residueid_c < residueid_d and residueid_d < residueid_e and residueid_e < residueid_f;


drop table EDH_site_2;
create table EDH_site_2 as 
	select f2.* from
	(select t2.* from
        (select a2.* from (select distinct a2.* from EDH_site_result_2 a2
        left join EDH_site_result_3 a3 on a2.pdbfileid=a3.pdbfileid and (a2.residueid_a=a3.residueid_a or a2.residueid_a=a3.residueid_b)
        where a3.pdbfileid is null) a2
        left join EDH_site_result_4 a4 on a2.pdbfileid=a4.pdbfileid and (a2.residueid_a=a4.residueid_a or a2.residueid_a=a4.residueid_b or a2.residueid_a=a4.residueid_c)
        where a4.pdbfileid is null) t2 
	left join EDH_site_result_5 a5 on t2.pdbfileid=a5.pdbfileid and (t2.residueid_a=a5.residueid_a or t2.residueid_a=a5.residueid_b or t2.residueid_a=a5.residueid_c or t2.residueid_a=a5.residueid_d)  
	where a5.pdbfileid is null) f2
	left join EDH_site_result_6 a6 on f2.pdbfileid=a6.pdbfileid and (f2.residueid_a=a6.residueid_a or f2.residueid_a=a6.residueid_b or f2.residueid_a=a6.residueid_c or f2.residueid_a=a6.residueid_d or f2.residueid_a=a6.residueid_e)
        where a6.pdbfileid is null;


drop table EDH_site_3;
create table EDH_site_3 as 
	select f3.* from
        (select t3.* from 
	(select a3.* from EDH_site_result_3 a3
        left join EDH_site_result_4 a4 on a3.pdbfileid=a4.pdbfileid and (a3.residueid_a=a4.residueid_a or a3.residueid_a=a4.residueid_b)
        where a4.pdbfileid is NULL) t3
	left join EDH_site_result_5 a5 on t3.pdbfileid=a5.pdbfileid and (t3.residueid_a=a5.residueid_a or t3.residueid_a=a5.residueid_b or t3.residueid_a=a5.residueid_c)
        where a5.pdbfileid is null) f3
	left join EDH_site_result_6 a6 on f3.pdbfileid=a6.pdbfileid and (f3.residueid_a=a6.residueid_a or f3.residueid_a=a6.residueid_b or f3.residueid_a=a6.residueid_c or f3.residueid_a=a6.residueid_d)
        where a6.pdbfileid is null;


drop table EDH_site_4;
create table EDH_site_4 as
        select t4.* from
        (select a4.* from EDH_site_result_4 a4
        left join EDH_site_result_5 a5 on a4.pdbfileid=a5.pdbfileid and (a4.residueid_a=a5.residueid_a or a4.residueid_a=a5.residueid_b)
        where a5.pdbfileid is NULL) t4
        left join EDH_site_result_6 a6 on t4.pdbfileid=a6.pdbfileid and (t4.residueid_a=a6.residueid_a or t4.residueid_a=a6.residueid_b or t4.residueid_a=a6.residueid_c)
        where a6.pdbfileid is null;

drop table EDH_site_5;
create table EDH_site_5 as
        select a5.* from EDH_site_result_5 a5
        left join EDH_site_result_6 a6 on a5.pdbfileid=a6.pdbfileid and (a5.residueid_a=a6.residueid_a or a5.residueid_a=a6.residueid_b)
        where a6.pdbfileid is NULL;


drop table EDH_site23456;
create table EDH_site23456 as 
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
                 chainid_d,resname_d,resseq_d,residueid_d,
		 chainid_e,resname_e,resseq_e,residueid_e,
		 chainid_f,resname_f,resseq_f,residueid_f,
                 atomid_a,atomid_b,atomid_c,atomid_d,atomid_e,atomid_f,
                atomname_a,atomname_b,atomname_c,atomname_d,atomname_e,atomname_f,
                contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd,contact_ae,contact_be,contact_ce,contact_de,contact_af,contact_bf,contact_cf,contact_df,contact_ef,6 as site_count
                from EDH_site_result_6
	union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
                 chainid_d,resname_d,resseq_d,residueid_d,
                 chainid_e,resname_e,resseq_e,residueid_e,
                'NaN' as chainid_f,'NaN' as resname_f,-9999 as resseq_f,-9999 as residueid_f,
                 atomid_a,atomid_b,atomid_c,atomid_d,atomid_e,-9999 as atomid_f,
                atomname_a,atomname_b,atomname_c,atomname_d,atomname_e,'NaN' as atomname_f,
		contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd,contact_ae,contact_be,contact_ce,contact_de,cast(1111 as bit) as contact_af,cast(1111 as bit) as contact_bf,cast(1111 as bit) as contact_cf,cast(1111 as bit) as contact_df,cast(1111 as bit) as contact_ef, 5 as site_count
                from EDH_site_5	
	union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
                 chainid_d,resname_d,resseq_d,residueid_d,
		 'NaN' as chainid_e,'NaN' as resname_e,-9999 as resseq_e,-9999 as residueid_e,
                 'NaN' as chainid_f,'NaN' as resname_f,-9999 as resseq_f,-9999 as residueid_f,
                 atomid_a,atomid_b,atomid_c,atomid_d,-9999 as atomid_e,-9999 as atomid_f,
                atomname_a,atomname_b,atomname_c,atomname_d,'NaN' as atomname_e,'NaN' as atomname_f,
		contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd,cast(1111 as bit) as contact_ae,cast(1111 as bit) as contact_be,cast(1111 as bit) as contact_ce,cast(1111 as bit) as contact_de,cast(1111 as bit) as contact_af,cast(1111 as bit) as contact_bf,cast(1111 as bit) as contact_cf,cast(1111 as bit) as contact_df,cast(1111 as bit) as contact_ef, 4 as site_count
                from EDH_site_4	
	union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 chainid_c,resname_c,resseq_c,residueid_c,
		 'NaN' as chainid_d,'NaN' as resname_d,-9999 as resseq_d,-9999 as residueid_d,
                 'NaN' as chainid_e,'NaN' as resname_e,-9999 as resseq_e,-9999 as residueid_e,
                 'NaN' as chainid_f,'NaN' as resname_f,-9999 as resseq_f,-9999 as residueid_f,
                 atomid_a,atomid_b,atomid_c,-9999 as atomid_d,-9999 as atomid_e,-9999 as atomid_f,
                atomname_a,atomname_b,atomname_c,'NaN' as atomname_d,'NaN' as atomname_e,'NaN' as atomname_f,
		contact_ab,contact_bc,contact_ca,cast(1111 as bit) as contact_ad,cast(1111 as bit) as contact_bd,cast(1111 as bit) as contact_cd,cast(1111 as bit) as contact_ae,cast(1111 as bit) as contact_be,cast(1111 as bit) as contact_ce,cast(1111 as bit) as contact_de,cast(1111 as bit) as contact_af,cast(1111 as bit) as contact_bf,cast(1111 as bit) as contact_cf,cast(1111 as bit) as contact_df,cast(1111 as bit) as contact_ef, 3 as site_count
                from EDH_site_3
	union
    select pdbfileid,chainid_a,resname_a,resseq_a,residueid_a,
                 chainid_b,resname_b,resseq_b,residueid_b,
                 'NaN' as chainid_c,'NaN' as resname_c,-9999 as resseq_c,-9999 as residueid_c,
                 'NaN' as chainid_d,'NaN' as resname_d,-9999 as resseq_d,-9999 as residueid_d,
                 'NaN' as chainid_e,'NaN' as resname_e,-9999 as resseq_e,-9999 as residueid_e,
                 'NaN' as chainid_f,'NaN' as resname_f,-9999 as resseq_f,-9999 as residueid_f,
                 atomid_a,atomid_b,-9999 as atomid_c,-9999 as atomid_d,-9999 as atomid_e,-9999 as atomid_f,
                atomname_a,atomname_b,'NaN' as atomname_c,'NaN' as atomname_d,'NaN' as atomname_e,'NaN' as atomname_f,
		contact_ab,cast(1111 as bit) as contact_bc,cast(1111 as bit) as contact_ca,cast(1111 as bit) as contact_ad,cast(1111 as bit) as contact_bd,cast(1111 as bit) as contact_cd,cast(1111 as bit) as contact_ae,cast(1111 as bit) as contact_be,cast(1111 as bit) as contact_ce,cast(1111 as bit) as contact_de,cast(1111 as bit) as contact_af,cast(1111 as bit) as contact_bf,cast(1111 as bit) as contact_cf,cast(1111 as bit) as contact_df,cast(1111 as bit) as contact_ef, 2 as site_count
                from EDH_site_2;


drop table metal_predict_sites;
create table metal_predict_sites as
        select t1.*, t2.x as x_a, t2.y as y_a, t2.z as z_a,
        t3.x as x_b, t3.y as y_b, t3.z as z_b,
        t4.x as x_c, t4.y as y_c, t4.z as z_c,
        t5.x as x_d, t5.y as y_d, t5.z as z_d,
	t5.x as x_e, t5.y as y_e, t5.z as z_e,
	t5.x as x_f, t5.y as y_f, t5.z as z_f
        from (select distinct b.pdbid,a.* from EDH_site23456 a
        left join neighborhood.pdbfiles b
        on a.pdbfileid=b.pdbfileid) t1
        left join neighborhood.atoms t2
        on t1.pdbfileid = t2.pdbfileid and t1.residueid_a = t2.residueid and t1.atomid_a = t2.atomid
         left join neighborhood.atoms t3
        on t1.pdbfileid = t3.pdbfileid and t1.residueid_b = t3.residueid and t1.atomid_b = t3.atomid
         left join neighborhood.atoms t4
        on t1.pdbfileid = t4.pdbfileid and t1.residueid_c = t4.residueid and t1.atomid_c = t4.atomid
         left join neighborhood.atoms t5
        on t1.pdbfileid = t5.pdbfileid and t1.residueid_d = t5.residueid and t1.atomid_d = t5.atomid
	 left join neighborhood.atoms t6
        on t1.pdbfileid = t6.pdbfileid and t1.residueid_e = t6.residueid and t1.atomid_e = t6.atomid
	 left join neighborhood.atoms t7
        on t1.pdbfileid = t7.pdbfileid and t1.residueid_f = t7.residueid and t1.atomid_f = t7.atomid;


----delete from metal_predict_sites where chainid_a != chainid_b and site_count=2;
delete from metal_predict_sites where atomname_a in ('-CG-') or atomname_b in ('-CG-')
      or atomname_c in ('-CG-') or atomname_d in ('-CG-') ;


alter table metal_predict_sites add column id integer;
update metal_predict_sites a set
        id = b.num from (SELECT ROW_NUMBER() over(ORDER bY pdbid,residueid_a,residueid_b,residueid_c,residueid_d,residueid_e,residueid_f DESC ) As num,
                pdbid,residueid_a,residueid_b,residueid_c,residueid_d,residueid_e,residueid_f from metal_predict_sites) b
        where a.pdbid=b.pdbid and a.residueid_a=b.residueid_a and a.residueid_b=b.residueid_b and a.residueid_c=b.residueid_c and a.residueid_d=b.residueid_d and a.residueid_e=b.residueid_e and a.residueid_f=b.residueid_f;

alter table metal_predict_sites add column resi_type text;

update metal_predict_sites set resi_type = 'E_D' where 
	(resname_a in ('GLU','ASP') and resname_b in ('GLU','ASP') and site_count=2) or 
	(resname_a in ('GLU','ASP') and resname_b in ('GLU','ASP') and resname_c in ('GLU','ASP') and site_count=3) or 
	(resname_a in ('GLU','ASP') and resname_b in ('GLU','ASP') and resname_c in ('GLU','ASP') and resname_d in ('GLU','ASP') and site_count=4) or
	(resname_a in ('GLU','ASP') and resname_b in ('GLU','ASP') and resname_c in ('GLU','ASP') and resname_d in ('GLU','ASP') and resname_e in ('GLU','ASP') and site_count=5) or
	(resname_a in ('GLU','ASP') and resname_b in ('GLU','ASP') and resname_c in ('GLU','ASP') and resname_d in ('GLU','ASP') and resname_e in ('GLU','ASP') and resname_f in ('GLU','ASP') and site_count=6);

update metal_predict_sites set resi_type = 'H_H' where
        (resname_a in ('HIS') and resname_b in ('HIS') and site_count=2) or
        (resname_a in ('HIS') and resname_b in ('HIS') and resname_c in ('HIS') and site_count=3) or
        (resname_a in ('HIS') and resname_b in ('HIS') and resname_c in ('HIS') and resname_d in ('HIS') and site_count=4) or
        (resname_a in ('HIS') and resname_b in ('HIS') and resname_c in ('HIS') and resname_d in ('HIS') and resname_e in ('HIS') and site_count=5) or
        (resname_a in ('HIS') and resname_b in ('HIS') and resname_c in ('HIS') and resname_d in ('HIS') and resname_e in ('HIS') and resname_f in ('HIS') and site_count=6);

update metal_predict_sites set resi_type = 'EDH' where resi_type is null;

------ 5. Add atomic coordinate column (EDH)


drop table pre_metal_coord_sites;
create table pre_metal_coord_sites as
        select distinct id,pdbid, pdbfileid, 
	resname_a, residueid_a,
	resname_b, residueid_b, 
	resname_c, residueid_c,
	resname_d, residueid_d,
	resname_e, residueid_e,
	resname_f, residueid_f,
	contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd,contact_ae,contact_be,contact_ce,contact_de,contact_af,contact_bf,contact_cf,contact_df,contact_ef,
        resi_type,site_count from metal_predict_sites; 


-----OD1 OD2 OE1 OE2 (a)
    --------ASP	
alter table pre_metal_coord_sites add column OD1_coords_a text;
update pre_metal_coord_sites t1 set OD1_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t1.resname_a='ASP';

alter table pre_metal_coord_sites add column OD2_coords_a text;
update pre_metal_coord_sites t1 set OD2_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t1.resname_a='ASP';

    ------GLU
update pre_metal_coord_sites t1 set OD1_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t1.resname_a='GLU';

update pre_metal_coord_sites t1 set OD2_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid and t1.resname_a='GLU';

-----OD1 OD2 OE1 OE2 (b)

alter table pre_metal_coord_sites add column OD1_coords_b text;
update pre_metal_coord_sites t1 set OD1_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t1.resname_b='ASP';

alter table pre_metal_coord_sites add column OD2_coords_b text;
update pre_metal_coord_sites t1 set OD2_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t1.resname_b='ASP';

    ------GLU
update pre_metal_coord_sites t1 set OD1_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t1.resname_b='GLU';

update pre_metal_coord_sites t1 set OD2_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid and t1.resname_b='GLU';

-----OD1 OD2 OE1 OE2 (c)
alter table pre_metal_coord_sites add column OD1_coords_c text;
update pre_metal_coord_sites t1 set OD1_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t1.resname_c='ASP';

alter table pre_metal_coord_sites add column OD2_coords_c text;
update pre_metal_coord_sites t1 set OD2_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t1.resname_c='ASP';

    ------GLU
update pre_metal_coord_sites t1 set OD1_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t1.resname_c='GLU';

update pre_metal_coord_sites t1 set OD2_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid and t1.resname_c='GLU';

-----OD1 OD2 OE1 OE2 (d)
alter table pre_metal_coord_sites add column OD1_coords_d text;
update pre_metal_coord_sites t1 set OD1_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t1.resname_d='ASP';

alter table pre_metal_coord_sites add column OD2_coords_d text;
update pre_metal_coord_sites t1 set OD2_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t1.resname_d='ASP';

    ------GLU
update pre_metal_coord_sites t1 set OD1_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t1.resname_d='GLU';

update pre_metal_coord_sites t1 set OD2_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid and t1.resname_d='GLU';

-----OD1 OD2 OE1 OE2 (e)
alter table pre_metal_coord_sites add column OD1_coords_e text;
update pre_metal_coord_sites t1 set OD1_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid and t1.resname_e='ASP';

alter table pre_metal_coord_sites add column OD2_coords_e text;
update pre_metal_coord_sites t1 set OD2_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid and t1.resname_e='ASP';

    ------GLU
update pre_metal_coord_sites t1 set OD1_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid and t1.resname_e='GLU';

update pre_metal_coord_sites t1 set OD2_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid and t1.resname_e='GLU';

-----OD1 OD2 OE1 OE2 (f)
alter table pre_metal_coord_sites add column OD1_coords_f text;
update pre_metal_coord_sites t1 set OD1_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid and t1.resname_f='ASP';

alter table pre_metal_coord_sites add column OD2_coords_f text;
update pre_metal_coord_sites t1 set OD2_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid and t1.resname_f='ASP';

    ------GLU
update pre_metal_coord_sites t1 set OD1_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid and t1.resname_f='GLU';

update pre_metal_coord_sites t1 set OD2_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-OE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid and t1.resname_f='GLU';


----- CD2 coord
alter table pre_metal_coord_sites add column CD2_coords_a text;
update pre_metal_coord_sites t1 set CD2_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column CD2_coords_b text;
update pre_metal_coord_sites t1 set CD2_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column CD2_coords_c text;
update pre_metal_coord_sites t1 set CD2_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column CD2_coords_d text;
update pre_metal_coord_sites t1 set CD2_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column CD2_coords_e text;
update pre_metal_coord_sites t1 set CD2_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column CD2_coords_f text;
update pre_metal_coord_sites t1 set CD2_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------ CE1 coord
alter table pre_metal_coord_sites add column CE1_coords_a text;
update pre_metal_coord_sites t1 set CE1_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column CE1_coords_b text;
update pre_metal_coord_sites t1 set CE1_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column CE1_coords_c text;
update pre_metal_coord_sites t1 set CE1_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column CE1_coords_d text;
update pre_metal_coord_sites t1 set CE1_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column CE1_coords_e text;
update pre_metal_coord_sites t1 set CE1_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column CE1_coords_f text;
update pre_metal_coord_sites t1 set CE1_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CE1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- N coord
alter table pre_metal_coord_sites add column N_coords_a text;
update pre_metal_coord_sites t1 set N_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column N_coords_b text;
update pre_metal_coord_sites t1 set N_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column N_coords_c text;
update pre_metal_coord_sites t1 set N_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column N_coords_d text;
update pre_metal_coord_sites t1 set N_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column N_coords_e text;
update pre_metal_coord_sites t1 set N_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column N_coords_f text;
update pre_metal_coord_sites t1 set N_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-N--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- ND1 coord
alter table pre_metal_coord_sites add column ND1_coords_a text;
update pre_metal_coord_sites t1 set ND1_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column ND1_coords_b text;
update pre_metal_coord_sites t1 set ND1_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column ND1_coords_c text;
update pre_metal_coord_sites t1 set ND1_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column ND1_coords_d text;
update pre_metal_coord_sites t1 set ND1_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column ND1_coords_e text;
update pre_metal_coord_sites t1 set ND1_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column ND1_coords_f text;
update pre_metal_coord_sites t1 set ND1_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-ND1') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- SG coord
alter table pre_metal_coord_sites add column SG_coords_a text;
update pre_metal_coord_sites t1 set SG_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column SG_coords_b text;
update pre_metal_coord_sites t1 set SG_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column SG_coords_c text;
update pre_metal_coord_sites t1 set SG_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column SG_coords_d text;
update pre_metal_coord_sites t1 set SG_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column SG_coords_e text;
update pre_metal_coord_sites t1 set SG_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column SG_coords_f text;
update pre_metal_coord_sites t1 set SG_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-SG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- O coord
alter table pre_metal_coord_sites add column O_coords_a text;
update pre_metal_coord_sites t1 set O_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column O_coords_b text;
update pre_metal_coord_sites t1 set O_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column O_coords_c text;
update pre_metal_coord_sites t1 set O_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column O_coords_d text;
update pre_metal_coord_sites t1 set O_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column O_coords_e text;
update pre_metal_coord_sites t1 set O_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column O_coords_f text;
update pre_metal_coord_sites t1 set O_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-O--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- NE2 coord
alter table pre_metal_coord_sites add column NE2_coords_a text;
update pre_metal_coord_sites t1 set NE2_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column NE2_coords_b text;
update pre_metal_coord_sites t1 set NE2_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column NE2_coords_c text;
update pre_metal_coord_sites t1 set NE2_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column NE2_coords_d text;
update pre_metal_coord_sites t1 set NE2_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column NE2_coords_e text;
update pre_metal_coord_sites t1 set NE2_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column NE2_coords_f text;
update pre_metal_coord_sites t1 set NE2_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-NE2') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

-------- C coord
alter table pre_metal_coord_sites add column C_coords_a text;
update pre_metal_coord_sites t1 set C_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column C_coords_b text;
update pre_metal_coord_sites t1 set C_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column C_coords_c text;
update pre_metal_coord_sites t1 set C_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column C_coords_d text;
update pre_metal_coord_sites t1 set C_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column C_coords_e text;
update pre_metal_coord_sites t1 set C_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column C_coords_f text;
update pre_metal_coord_sites t1 set C_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-C--') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- CA coord
alter table pre_metal_coord_sites add column CA_coords_a text;
update pre_metal_coord_sites t1 set CA_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column CA_coords_b text;
update pre_metal_coord_sites t1 set CA_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column CA_coords_c text;
update pre_metal_coord_sites t1 set CA_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column CA_coords_d text;
update pre_metal_coord_sites t1 set CA_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column CA_coords_e text;
update pre_metal_coord_sites t1 set CA_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column CA_coords_f text;
update pre_metal_coord_sites t1 set CA_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CA-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- CB coord
alter table pre_metal_coord_sites add column CB_coords_a text;
update pre_metal_coord_sites t1 set CB_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column CB_coords_b text;
update pre_metal_coord_sites t1 set CB_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column CB_coords_c text;
update pre_metal_coord_sites t1 set CB_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column CB_coords_d text;
update pre_metal_coord_sites t1 set CB_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column CB_coords_e text;
update pre_metal_coord_sites t1 set CB_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column CB_coords_f text;
update pre_metal_coord_sites t1 set CB_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CB-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- CG coord
alter table pre_metal_coord_sites add column CG_coords_a text;
update pre_metal_coord_sites t1 set CG_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column CG_coords_b text;
update pre_metal_coord_sites t1 set CG_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column CG_coords_c text;
update pre_metal_coord_sites t1 set CG_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column CG_coords_d text;
update pre_metal_coord_sites t1 set CG_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column CG_coords_e text;
update pre_metal_coord_sites t1 set CG_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column CG_coords_f text;
update pre_metal_coord_sites t1 set CG_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CG-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;

------- CD coord
alter table pre_metal_coord_sites add column CD_coords_a text;
update pre_metal_coord_sites t1 set CD_coords_a=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_a=t2.residueid;

alter table pre_metal_coord_sites add column CD_coords_b text;
update pre_metal_coord_sites t1 set CD_coords_b=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_b=t2.residueid;

alter table pre_metal_coord_sites add column CD_coords_c text;
update pre_metal_coord_sites t1 set CD_coords_c=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_c=t2.residueid;

alter table pre_metal_coord_sites add column CD_coords_d text;
update pre_metal_coord_sites t1 set CD_coords_d=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_d=t2.residueid;

alter table pre_metal_coord_sites add column CD_coords_e text;
update pre_metal_coord_sites t1 set CD_coords_e=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_e=t2.residueid;

alter table pre_metal_coord_sites add column CD_coords_f text;
update pre_metal_coord_sites t1 set CD_coords_f=t2.atom_coords
        from (select pdbfileid,residueid,atomid,atomname, concat(x,',',y,',',z) as atom_coords from neighborhood.atoms where atomname = '-CD-') t2
        where t1.pdbfileid=t2.pdbfileid and t1.residueid_f=t2.residueid;
