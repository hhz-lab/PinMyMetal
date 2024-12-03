---- get chainid_ion
drop table pre_chain;
create table pre_chain (
        sitetype char(3),
        id integer,
        chainid char(3));
\copy pre_chain from './get_premetal_chainid.csv'

alter table metal_predict_sites_2 drop column chainid_ion;
alter table metal_predict_sites_2 add column chainid_ion char(3);
update metal_predict_sites_2 t1 set chainid_ion=t2.chainid from pre_chain t2 where t1.id=t2.id and t2.sitetype='edh';

alter table zncu_predict_sites_2 drop column chainid_ion;
alter table zncu_predict_sites_2 add column chainid_ion char(3);
update zncu_predict_sites_2 t1 set chainid_ion=t2.chainid from pre_chain t2 where t1.id=t2.id and t2.sitetype='ch';

--- Find the nearest prediction site for each experimental metal site in the prediction dataset
drop table match_exp2pre;
create table match_exp2pre (
        exp_sitetype char(3),
        exp_pdbid char(4),
        exp_id integer,
        exp_metal_x float,
        exp_metal_y float,
        exp_metal_z float,
        exp_metalname char(3),
        pred_sitetype char(3),
        pred_pdbid char(4),
        pred_id integer,
        pred_metal_x float,
        pred_metal_y float,
        pred_metal_z float,
        pred_proba_ismetal float,
        distance float);
\copy match_exp2pre from './match_exp2pre.csv' delimiter E'\t' csv header

alter table match_exp2pre add column pre_is_exp bool default false;
update match_exp2pre t1 set pre_is_exp=true from
        (select distinct on (exp_pdbid,exp_id) exp_pdbid,exp_id, distance from
        match_exp2pre order by exp_pdbid,exp_id, distance) t2
        where t1.exp_pdbid=t2.exp_pdbid and t1.exp_id=t2.exp_id and t1.distance=t2.distance and t1.distance < 2.5;


------ The CH sites, EDH sites and experimental sites were clustered according to dist < 2.5
drop table metals_cluster;
create table metals_cluster
        (pdbid char(4),
        clusterid real,
        id real,
        sitetype char(3),
        ligands text,
        is_rep bool,
        atomname_ion text,
        exp_site bool,
        dist_pre_exp float);
\copy metals_cluster from 'metalsites_cluster.csv' delimiter E'\t' csv header


-----pre is exp sites
drop table pre_sites;
create table pre_sites as 
	select distinct a.pdbid,metal_x,metal_y,metal_z, chainid_ion, a.pre_id, a.sitetype,
	chainid_a, chainid_b,chainid_c,chainid_d,'NaN' as chainid_e, 'NaN' as chainid_f, 
	resseq_a, resseq_b,resseq_c,resseq_d, -9999 as resseq_e, -9999 as resseq_f, residueid_a, residueid_b,
	atom_a,atom_b,atom_c,atom_d,'NaN' as atom_e, 'NaN' as atom_f,
	resname_a,resname_b,resname_c,resname_d, 'NaN' as resname_e, 'NaN' as resname_f,
	dist_am,dist_bm,dist_cm,dist_dm, 0 as dist_em, 0 as dist_fm, proba_ismetal,site_count
	from (select distinct pdbid,id as pre_id,sitetype from metals_cluster where is_rep is true and sitetype='ch') a
    left join zncu_predict_sites_2 b on a.pre_id=b.id
    union
	select distinct c.pdbid,metal_x,metal_y,metal_z, chainid_ion, c.pre_id, c.sitetype,
       chainid_a, chainid_b,chainid_c,chainid_d,chainid_e,chainid_f,
        resseq_a, resseq_b,resseq_c,resseq_d,resseq_e,resseq_f, residueid_a, residueid_b,
        atom_a,atom_b,atom_c,atom_d,atom_e,atom_f,
	resname_a,resname_b,resname_c,resname_d,resname_e,resname_f,
        dist_am,dist_bm,dist_cm,dist_dm,dist_em,dist_fm, proba_ismetal,site_count
        from (select distinct pdbid,id as pre_id,sitetype from metals_cluster where is_rep is true and sitetype='edh') c
    left join metal_predict_sites_2 d on c.pre_id=d.id;

alter table pre_sites add column new_preid integer;
update pre_sites a set
        new_preid = b.num from (SELECT ROW_NUMBER() over(ORDER bY pre_id,sitetype DESC ) As num,
                pre_id,sitetype from pre_sites) b
        where a.pre_id=b.pre_id and a.sitetype=b.sitetype;

alter table match_exp2pre add column new_preid integer;
update match_exp2pre t1 set
	new_preid = t2.new_preid from (select distinct pre_id,sitetype,new_preid from pre_sites) t2
	where t1.pred_id=t2.pre_id and t1.pred_sitetype=t2.sitetype;

----- all exp sites
drop table exp_sites;
create table exp_sites as select a.*, chainid,resseq,chainids_lig,resseqs_lig,resnames_lig
	from (select distinct id as exp_id, pdbid, residueid_ion,x_ion,y_ion,z_ion,bench,resname_ion,which_metal from transition_metal_coord
		where resname_ion in ('_MN','FE2','_FE','_CO','_NI','_CU','_ZN')) a
	left join all_metal b
	on a.pdbid=b.pdbid and a.residueid_ion=b.residueid_ion;

alter table exp_sites add column new_preid smallint default -1;
update exp_sites t1 set new_preid=t2.new_preid from
	(select * from match_exp2pre where pre_is_exp is true) t2
	where t1.exp_id=t2.exp_id;

------ Labeling and classifying predicted sites
	---- add bench column
alter table pre_sites add column bench text;
update pre_sites t1 set bench = 'true' from
       (select distinct exp_id,new_preid,bench from exp_sites where bench is true) t2
                 where t1.new_preid=t2.new_preid;
update pre_sites t1 set bench = 'false' from
       (select distinct exp_id,new_preid,bench from exp_sites where bench is false) t2
                 where t1.new_preid=t2.new_preid;
update pre_sites t1 set bench = 'predict' where bench is null;

---- The variable exp_tag=1 indicates that the predicted site is a true experimental site from the benchmark.
alter table pre_sites add column exp_tag smallint default 0;
update pre_sites set exp_tag = 1 where bench = 'true';

alter table exp_sites add column metal_label text;
update exp_sites set metal_label = 'Manganese' where resname_ion = '_MN';
update exp_sites set metal_label = 'Iron' where resname_ion in ('FE2','_FE');
update exp_sites set metal_label = 'Cobalt' where resname_ion = '_CO';
update exp_sites set metal_label = 'Nickel' where resname_ion = '_NI';
update exp_sites set metal_label = 'Copper' where resname_ion = '_CU';
update exp_sites set metal_label = 'Zinc' where resname_ion = '_ZN';

