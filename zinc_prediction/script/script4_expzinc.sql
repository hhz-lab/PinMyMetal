drop table hydrophobic_proba;
create table hydrophobic_proba
        (id integer,
        pearson_c float,
        pearson_s float);
\copy hydrophobic_proba from './hydro_proba.csv' delimiter ','

alter table hydrophobic_pre drop column zinc_tag;
alter table hydrophobic_pre add column zinc_tag bool default false;
update hydrophobic_pre t1 set zinc_tag=true from (select distinct id, (pearson_c+pearson_s)/2 as p_pearson from hydrophobic_proba) t2
        where t1.id=t2.id and p_pearson > 0.5;

alter table zinc_predict_site234_2 add zinc_tag bool default false;
update zinc_predict_site234_2 t1 set zinc_tag=t2.zinc_tag from hydrophobic_pre t2 where t1.id=t2.id;
update zinc_predict_site234_2 t1 set selected_site is false where zinc_tag is false and bench is false and site_count >=3;

alter table zinc_predict_site234_2 add p_value float;
update zinc_predict_site234_2 t1 set p_value=t2.p_pearson from (select distinct id, (pearson_c+pearson_s)/2 as p_pearson from hydrophobic_proba) t2
        where t1.id=t2.id;

drop table ml_predict_result;
create table ml_predict_result 
	(id int,
	result smallint);
\copy ml_predict_result from './ml_result.csv' delimiter ',' csv header

alter table zinc_predict_site234_2 add column ml_result smallint;

update zinc_predict_site234_2 t1 set ml_result=t2.result from ml_predict_result t2 where t1.id=t2.id;

update zinc_predict_site234_2 t1 set ml_result=1 where zinc_tag is true and site_count >=3;

drop table mvc_proba;
create table mvc_proba
        (id int,
        proba float);
\copy mvc_proba from './mvc_proba.csv' delimiter ',' csv header

update zinc_predict_site234_2 t1 set p_value=t2.proba from mvc_proba t2 where t1.id=t2.id;

drop table exp_pre_site;
create table exp_pre_site as
        select distinct pdbid,'pre' as site_type, conc_comma(id||'_'||zinc_x||'_'||zinc_y||'_'||zinc_z) as zincs  from
    (select distinct pdbid,id,zinc_x,zinc_y,zinc_z from zinc_predict_site234_2) a group by pdbid
        union select distinct pdbid, 'exp' as site_type, conc_comma(residueid_ion||'_'||x_ion||'_'||y_ion||'_'||z_ion) as zincs from
        (select distinct pdbid,residueid_ion,x_ion,y_ion,z_ion from pdb_zinc_coordinate) a group by pdbid order by pdbid ;

alter table exp_pre_site add column redundancy bool default false;
update exp_pre_site t1 set redundancy = true from (select pdbid,count(*) as pdbid_count from exp_pre_site group by pdbid ) t2 where t1.pdbid=t2.pdbid and t2.pdbid_count=2;

drop table exp_zinc_site2;
create table exp_zinc_site2 as
        select a.pdbid,a.residueid_ion,a.x_ion,a.y_ion,a.z_ion,
        a.residueid_lig as residueid_a, b.residueid_lig as residueid_b,a.bench from
        (select * from pdb_zinc_coordinate where site_count = 2 and resname_lig in ('CYS','HIS') ) a
        left join (select * from pdb_zinc_coordinate where site_count = 2 and resname_lig in ('CYS','HIS')) b
        on a.pdbfileid = b.pdbfileid and a.residueid_ion=b.residueid_ion
        where a.residueid_lig < b.residueid_lig and a.atomid_lig < b.atomid_lig;

drop table exp_zinc_site3;
create table exp_zinc_site3 as
        select a.pdbid,a.residueid_ion,a.x_ion,a.y_ion,a.z_ion,
        a.residueid_lig as residueid_a, b.residueid_lig as residueid_b,
        c.residueid_lig as residueid_c,a.bench  from
        (select * from pdb_zinc_coordinate where site_count = 3 and resname_lig in ('CYS','HIS')) a
        left join (select * from pdb_zinc_coordinate where site_count = 3 and resname_lig in ('CYS','HIS')) b
        on a.pdbfileid = b.pdbfileid and a.residueid_ion=b.residueid_ion
        left join (select * from pdb_zinc_coordinate where site_count = 3 and resname_lig in ('CYS','HIS')) c
        on a.pdbfileid = c.pdbfileid and a.residueid_ion=c.residueid_ion
        where a.residueid_lig < b.residueid_lig and a.atomid_lig < b.atomid_lig and b.residueid_lig < c.residueid_lig and b.atomid_lig < c.atomid_lig;

drop table exp_zinc_site4;
create table exp_zinc_site4 as
        select a.pdbid,a.residueid_ion,a.x_ion,a.y_ion,a.z_ion,
        a.residueid_lig as residueid_a, b.residueid_lig as residueid_b,
        c.residueid_lig as residueid_c,
        d.residueid_lig as residueid_d,a.bench from
        (select * from pdb_zinc_coordinate where site_count = 4 and resname_lig in ('CYS','HIS')) a
        left join (select * from pdb_zinc_coordinate where site_count = 4 and resname_lig in ('CYS','HIS')) b
        on a.pdbfileid = b.pdbfileid and a.residueid_ion=b.residueid_ion
        left join (select * from pdb_zinc_coordinate where site_count = 4 and resname_lig in ('CYS','HIS')) c
        on a.pdbfileid = c.pdbfileid and a.residueid_ion=c.residueid_ion
        left join (select * from pdb_zinc_coordinate where site_count = 4 and resname_lig in ('CYS','HIS')) d
        on a.pdbfileid = d.pdbfileid and a.residueid_ion=d.residueid_ion
        where a.residueid_lig < b.residueid_lig and a.atomid_lig < b.atomid_lig and b.residueid_lig < c.residueid_lig and b.atomid_lig < c.atomid_lig
        and c.residueid_lig < d.residueid_lig and c.atomid_lig < d.atomid_lig;

drop table exp_site234;
create table exp_site234 as
        select pdbid, residueid_ion,x_ion,y_ion,z_ion,bench from exp_zinc_site2 union
        select pdbid, residueid_ion,x_ion,y_ion,z_ion,bench from exp_zinc_site3 union
        select pdbid, residueid_ion,x_ion,y_ion,z_ion,bench from exp_zinc_site4;

