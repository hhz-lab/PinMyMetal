drop table hydro_pre_chedh;
create table hydro_pre_chedh
        (id integer,
        tag float,
        c_value float,
        solv float,
        sitetype char(3));
\copy hydro_pre_chedh from './hydro_pre_chedh.csv' delimiter ','

drop table chemtype_pre_chedh;
create table chemtype_pre_chedh
        (id integer,
        chem_class float,
        distance float,
        sitetype char(3));
\copy chemtype_pre_chedh from './chem_pre_chedh.csv' delimiter ','

alter table hydro_pre_chedh add column resi_type text;
update hydro_pre_chedh t1 set resi_type = t2.resi_type from zncu_predict_sites_2 t2 where t1.id=t2.id and t1.sitetype='ch';
update hydro_pre_chedh t1 set resi_type = t2.resi_type from metal_predict_sites_2 t2 where t1.id=t2.id and t1.sitetype='edh';

alter table hydro_pre_chedh add column site_count smallint;
update hydro_pre_chedh t1 set site_count=t2.site_count from zncu_predict_sites_2 t2 where t1.id=t2.id and t1.sitetype='ch';
update hydro_pre_chedh t1 set site_count=t2.site_count from metal_predict_sites_2 t2 where t1.id=t2.id and t1.sitetype='edh';


alter table hydro_pre_chedh add column c_value_exp float;
update hydro_pre_chedh t1 set c_value_exp = t2.c_value from c_value_exp_chedh t2
        where t1.tag=t2.tag and t1.resi_type=t2.resi_type and t1.sitetype=t2.sitetype;

alter table hydro_pre_chedh add column solv_exp float;
update hydro_pre_chedh t1 set solv_exp = t2.solv from c_value_exp_chedh t2
        where t1.tag=t2.tag and t1.resi_type=t2.resi_type and t1.sitetype=t2.sitetype;

