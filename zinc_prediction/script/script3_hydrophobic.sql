drop table hydrophobic_pre;
create table hydrophobic_pre
        (id smallint,
        tag float,
        c_value float,
	solv float);
\copy hydrophobic_pre from './zinc_hydro.csv' delimiter ','

--alter table hydrophobic_pre add column atom_type text;
--update hydrophobic_pre t1 set atom_type=t2.atom_type from zinc_predict_site234_2 t2 where t1.id=t2.id;

alter table hydrophobic_pre add column resi_type text;
update hydrophobic_pre t1 set resi_type=t2.resi_type from zinc_predict_site234_2 t2 where t1.id=t2.id;

alter table hydrophobic_pre add column site_count smallint;
update hydrophobic_pre t1 set site_count=t2.site_count from zinc_predict_site234_2 t2 where t1.id=t2.id;

alter table hydrophobic_pre add column c_value_exp float;
update hydrophobic_pre t1 set c_value_exp = t2.c_value from c_value_exp t2
        where t1.tag=t2.tag and t1.resi_type=t2.resi_type and t1.site_count=t2.site_count;

alter table hydrophobic_pre add column solv_exp float;
update hydrophobic_pre t1 set solv_exp = t2.solv from c_value_exp t2
        where t1.tag=t2.tag and t1.resi_type=t2.resi_type and t1.site_count=t2.site_count;


drop table ml_pre;
create table ml_pre as
select * from (
select a.*,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21 from
        (select id, resname_a, resname_b,atomname_a,atomname_b,
        dist_ab as atom_dist,ca_dist,cb_dist,angle,angle1,angle2,angle3,atom_type
        from zinc_predict_site234_2 where site_count=2) a
        left join (select id, tag, solv as s1, c_value as c1 from hydrophobic_pre where tag=2) t1 on a.id=t1.id
        left join (select id, tag, solv as s2, c_value as c2 from hydrophobic_pre where tag=2.25) t2 on a.id=t2.id
        left join (select id, tag, solv as s3, c_value as c3 from hydrophobic_pre where tag=2.5) t3 on a.id=t3.id
        left join (select id, tag, solv as s4, c_value as c4 from hydrophobic_pre where tag=2.75) t4 on a.id=t4.id
        left join (select id, tag, solv as s5, c_value as c5 from hydrophobic_pre where tag=3) t5 on a.id=t5.id
        left join (select id, tag, solv as s6, c_value as c6 from hydrophobic_pre where tag=3.25) t6 on a.id=t6.id
        left join (select id, tag, solv as s7, c_value as c7 from hydrophobic_pre where tag=3.5) t7 on a.id=t7.id
        left join (select id, tag, solv as s8, c_value as c8 from hydrophobic_pre where tag=3.75) t8 on a.id=t8.id
        left join (select id, tag, solv as s9, c_value as c9 from hydrophobic_pre where tag=4) t9 on a.id=t9.id
        left join (select id, tag, solv as s10, c_value as c10 from hydrophobic_pre where tag=4.25) t10 on a.id=t10.id
        left join (select id, tag, solv as s11, c_value as c11 from hydrophobic_pre where tag=4.5) t11 on a.id=t11.id
        left join (select id, tag, solv as s12, c_value as c12 from hydrophobic_pre where tag=4.75) t12 on a.id=t12.id
        left join (select id, tag, solv as s13, c_value as c13 from hydrophobic_pre where tag=5) t13 on a.id=t13.id
        left join (select id, tag, solv as s14, c_value as c14 from hydrophobic_pre where tag=5.25) t14 on a.id=t14.id
        left join (select id, tag, solv as s15, c_value as c15 from hydrophobic_pre where tag=5.5) t15 on a.id=t15.id
        left join (select id, tag, solv as s16, c_value as c16 from hydrophobic_pre where tag=5.75) t16 on a.id=t16.id
        left join (select id, tag, solv as s17, c_value as c17 from hydrophobic_pre where tag=6) t17 on a.id=t17.id
        left join (select id, tag, solv as s18, c_value as c18 from hydrophobic_pre where tag=6.25) t18 on a.id=t18.id
        left join (select id, tag, solv as s19, c_value as c19 from hydrophobic_pre where tag=6.5) t19 on a.id=t19.id
        left join (select id, tag, solv as s20, c_value as c20 from hydrophobic_pre where tag=6.75) t20 on a.id=t20.id
        left join (select id, tag, solv as s21, c_value as c21 from hydrophobic_pre where tag=7) t21 on a.id=t21.id) a where c21 is not null;

alter table ml_pre add column a_ne2 smallint default 0;
alter table ml_pre add column a_nd1 smallint default 0;
alter table ml_pre add column b_ne2 smallint default 0;
alter table ml_pre add column b_nd1 smallint default 0;
alter table ml_pre add column a_ce1 smallint default 0;
alter table ml_pre add column a_cd2 smallint default 0;
alter table ml_pre add column b_ce1 smallint default 0;
alter table ml_pre add column b_cd2 smallint default 0;
alter table ml_pre add column a_sg smallint default 0;
alter table ml_pre add column b_sg smallint default 0;

update ml_pre set a_ne2 = 1 where atomname_a = '-NE2';
update ml_pre set a_nd1 = 1 where atomname_a = '-ND1';
update ml_pre set a_sg = 1 where atomname_a = '-SG-';
update ml_pre set a_ce1 = 1 where atomname_a = '-CE1';
update ml_pre set a_cd2 = 1 where atomname_a = '-CD2';

update ml_pre set b_ne2 = 1 where atomname_b = '-NE2';
update ml_pre set b_nd1 = 1 where atomname_b = '-ND1';
update ml_pre set b_sg = 1 where atomname_b = '-SG-';
update ml_pre set b_ce1 = 1 where atomname_b = '-CE1';
update ml_pre set b_cd2 = 1 where atomname_b = '-CD2';

update ml_pre set resname_a = 0 where resname_a = 'CYS';
update ml_pre set resname_a = 1 where resname_a = 'HIS';
update ml_pre set resname_b = 0 where resname_b = 'CYS';
update ml_pre set resname_b = 1 where resname_b = 'HIS';

