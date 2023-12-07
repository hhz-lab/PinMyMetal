drop table exp_pre_dist;
create table exp_pre_dist(
        pdbid char(4),
        exp_residueid integer,
        pre_id integer,
        dist float);
\copy exp_pre_dist from './exp_pre_dist.csv' delimiter ','

alter table exp_pre_dist add column bench text;
update exp_pre_dist t1 set bench=t2.bench from 
	(select distinct pdbid,residueid_ion,bench from pdb_zinc_coordinate) t2 
	where t1.pdbid=t2.pdbid and t1.exp_residueid=t2.residueid_ion;

alter table zinc_predict_site234_2 drop column selected_site;
alter table zinc_predict_site234_2 add column selected_site bool default false;
update zinc_predict_site234_2 t1 set selected_site=true from
        (select distinct on (pdbid,pre_id) pdbid,pre_id,exp_residueid,dist from (select * from exp_pre_dist where dist <= 2.5) a
                order by pdbid,pre_id,dist) t2
         where t1.pdbid=t2.pdbid and t1.id=t2.pre_id;

alter table zinc_predict_site234_2 add column bench text;
update zinc_predict_site234_2 t1 set bench=t2.bench from
        (select * from exp_pre_dist where pre_id in (select distinct id from zinc_predict_site234_2 where selected_site is true) )t2 where t1.id=t2.pre_id;


alter table exp_site234 add column id smallint default -1;
update exp_site234 t1 set id=t2.pre_id from 
	(select * from exp_pre_dist where pre_id in (select distinct id from zinc_predict_site234_2 where selected_site is true) ) t2 
	where t1.pdbid=t2.pdbid and t1.residueid_ion=t2.exp_residueid;

alter table zinc_predict_site234_2 add column exp_tag smallint default 0;
update zinc_predict_site234_2 set exp_tag = 1 where selected_site is true and bench = 'true';
