drop table hydro_proba_chedh;
create table hydro_proba_chedh
        (id integer,
        sitetype char(3),
        pearson_c float,
        pearson_s float,
        p_value float);
\copy hydro_proba_chedh from './hydro_proba.csv' delimiter ','
     ----- add metal_tag
alter table hydro_pre_chedh add column metal_tag bool default false;
update hydro_pre_chedh t1 set metal_tag=true from (select distinct id, sitetype, p_value from hydro_proba_chedh) t2
        where t1.id=t2.id and t1.sitetype=t2.sitetype and p_value > 0.5 and p_value !='NaN';

alter table metal_predict_sites_2 add metal_tag bool default false;
update metal_predict_sites_2 t1 set metal_tag=t2.metal_tag from hydro_pre_chedh t2 where t1.id=t2.id and t2.site_count >=4 and t2.sitetype='edh';

alter table zncu_predict_sites_2 add metal_tag bool default false;
update zncu_predict_sites_2 t1 set metal_tag=t2.metal_tag from hydro_pre_chedh t2 where t1.id=t2.id and t2.site_count >=3 and t2.sitetype='ch';

