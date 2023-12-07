----- for the pre_pre dist: script python3 pre_pre_dist.py > pre_pre_dist.csv
drop table pre_pre_dist;
create table pre_pre_dist(
        site_type text,
        pdbid char(4),
        id_a integer,
        id_b integer,
        dist float);
\copy pre_pre_dist from './pre_pre_dist.csv' delimiter ','

update zinc_predict_site234_2 set selected_site=true where id in (
select distinct id from (select distinct id_a as id from pre_pre_dist  where dist >2.5
        union select distinct id_b as id from pre_pre_dist  where dist >2.5) t
        where id not in (select distinct id_a as id from pre_pre_dist where dist <= 2.5 union select distinct id_b as id from pre_pre_dist where dist <= 2.5) );


