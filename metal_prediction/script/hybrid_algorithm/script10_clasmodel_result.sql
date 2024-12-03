drop table chedhclass_result;
create table chedhclass_result
        ( id int,
        sitetype char(3),
        metal_label smallint,
        proba_class float);
\copy chedhclass_result from 'chedh_metalclass_result.csv' delimiter ',' csv header

alter table pre_sites add column metal_label text;
update pre_sites t1 set metal_label='Manganese' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=1;
update pre_sites t1 set metal_label='FECONI' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=2;
update pre_sites t1 set metal_label='Copper' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=6;
update pre_sites t1 set metal_label='Zinc' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=7;

alter table pre_sites add column metal_pdb text;
update pre_sites t1 set metal_pdb='MN' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=1;
update pre_sites t1 set metal_pdb='8B' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=2;
update pre_sites t1 set metal_pdb='CU' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=6;
update pre_sites t1 set metal_pdb='ZN' from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype and t2.metal_label=7;


alter table pre_sites add column proba_class float;
update pre_sites t1 set proba_class=t2.proba_class from chedhclass_result t2 where t1.pre_id=t2.id and t1.sitetype=t2.sitetype;

