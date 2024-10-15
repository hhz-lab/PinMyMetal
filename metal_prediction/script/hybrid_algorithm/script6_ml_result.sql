drop table chmodel_result;
create table chmodel_result
        ( id int,
        result smallint,
        proba_ismetal float);
\copy chmodel_result from 'chmodel_result.csv' delimiter ',' csv header
--alter table zncu_predict_sites_2 drop column proba_ismetal;
alter table zncu_predict_sites_2 add column proba_ismetal float;
update zncu_predict_sites_2 t1 set proba_ismetal=t2.proba_ismetal from chmodel_result t2 where t1.id=t2.id;
update zncu_predict_sites_2 t1 set proba_ismetal=t2.p_value from hydro_proba_chedh t2 where t1.id=t2.id and t1.site_count >=3 and t2.sitetype='ch';

drop table edhmodel_result;
create table edhmodel_result
        ( id int,
        result smallint,
        proba_ismetal float);
\copy edhmodel_result from 'edhmodel_result.csv' delimiter ',' csv header
--alter table metal_predict_sites_2 drop column proba_ismetal;
alter table metal_predict_sites_2 add column proba_ismetal float;
update metal_predict_sites_2 t1 set proba_ismetal=t2.proba_ismetal from edhmodel_result t2 where t1.id=t2.id;
update metal_predict_sites_2 t1 set proba_ismetal=t2.p_value from hydro_proba_chedh t2 where t1.id=t2.id and t1.site_count >=4 and t2.sitetype='edh';


--alter table zncu_predict_sites_2 drop column mvc_result;
alter table zncu_predict_sites_2 add column mvc_result smallint default 0;
update zncu_predict_sites_2 set mvc_result=1 where proba_ismetal > 0.75;

--alter table metal_predict_sites_2 drop column mvc_result;
alter table metal_predict_sites_2 add column mvc_result smallint default 0;
update metal_predict_sites_2 t1 set mvc_result=1 where proba_ismetal > 0.75;

