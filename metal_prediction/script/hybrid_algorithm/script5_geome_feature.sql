drop table chedh_angel_dist;
create table chedh_angel_dist
        (id integer, sitetype char(3), ab_angle1 real, ab_angle2 real, ab_angle3 real, ac_angle1 real, ac_angle2 real, ac_angle3 real,
                ad_angle1 real, ad_angle2 real, ad_angle3 real,
            bc_angle1 real, bc_angle2 real, bc_angle3 real, bd_angle1 real, bd_angle2 real, bd_angle3 real,
        cd_angle1 real, cd_angle2 real, cd_angle3 real,
            ma_angle real, mb_angle real, mc_angle real, md_angle real,
            mo_angle_a real, mo_angle_b real,mo_angle_c real,mo_angle_d real,
            ca_dist_ab real, ca_dist_ac real, ca_dist_bc real,ca_dist_ad real, ca_dist_bd real, ca_dist_cd real,
            cb_dist_ab real,cb_dist_ac real, cb_dist_bc real, cb_dist_ad real, cb_dist_bd real, cb_dist_cd real,
            mCA_dist_a real, mCA_dist_b real, mCA_dist_c real, mCA_dist_d real,
            mCB_dist_a real, mCB_dist_b real, mCB_dist_c real, mCB_dist_d real,
            mo_dist_a real, mo_dist_b real, mo_dist_c real, mo_dist_d real);

\copy chedh_angel_dist from './calculate_angle_dist_site234.csv'

---- ch + edh features
drop table chedh_features;
create table chedh_features as
        select e.*, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21 from
        (select site_count, resname_a,resname_b,resname_c,resname_d,resi_type, b.*
        from (select distinct id, site_count, resname_a,resname_b,resname_c,resname_d,resi_type from zncu_predict_sites_2) a
        left join chedh_angel_dist b on a.id=b.id and b.sitetype='ch' union
        select site_count, resname_a,resname_b,resname_c,resname_d,resi_type, d.*
        from (select distinct id, site_count, resname_a,resname_b,resname_c,resname_d,resi_type from metal_predict_sites_2) c
        left join chedh_angel_dist d on c.id=d.id and d.sitetype='edh') e
        left join (select id, sitetype, tag, solv as s1, c_value as c1 from hydro_pre_chedh where tag=2) t1 on e.id=t1.id and e.sitetype=t1.sitetype
        left join (select id, sitetype, tag, solv as s2, c_value as c2 from hydro_pre_chedh where tag=2.25) t2 on e.id=t2.id and e.sitetype=t2.sitetype
        left join (select id, sitetype, tag, solv as s3, c_value as c3 from hydro_pre_chedh where tag=2.5) t3 on e.id=t3.id and e.sitetype=t3.sitetype
        left join (select id, sitetype, tag, solv as s4, c_value as c4 from hydro_pre_chedh where tag=2.75) t4 on e.id=t4.id and e.sitetype=t4.sitetype
        left join (select id, sitetype, tag, solv as s5, c_value as c5 from hydro_pre_chedh where tag=3) t5 on e.id=t5.id and e.sitetype=t5.sitetype
        left join (select id, sitetype, tag, solv as s6, c_value as c6 from hydro_pre_chedh where tag=3.25) t6 on e.id=t6.id and e.sitetype=t6.sitetype
        left join (select id, sitetype, tag, solv as s7, c_value as c7 from hydro_pre_chedh where tag=3.5) t7 on e.id=t7.id and e.sitetype=t7.sitetype
        left join (select id, sitetype, tag, solv as s8, c_value as c8 from hydro_pre_chedh where tag=3.75) t8 on e.id=t8.id and e.sitetype=t8.sitetype
        left join (select id, sitetype, tag, solv as s9, c_value as c9 from hydro_pre_chedh where tag=4) t9 on e.id=t9.id and e.sitetype=t9.sitetype
        left join (select id, sitetype, tag, solv as s10, c_value as c10 from hydro_pre_chedh where tag=4.25) t10 on e.id=t10.id and e.sitetype=t10.sitetype
        left join (select id, sitetype, tag, solv as s11, c_value as c11 from hydro_pre_chedh where tag=4.5) t11 on e.id=t11.id and e.sitetype=t11.sitetype
        left join (select id, sitetype, tag, solv as s12, c_value as c12 from hydro_pre_chedh where tag=4.75) t12 on e.id=t12.id and e.sitetype=t12.sitetype
        left join (select id, sitetype, tag, solv as s13, c_value as c13 from hydro_pre_chedh where tag=5) t13 on e.id=t13.id and e.sitetype=t13.sitetype
        left join (select id, sitetype, tag, solv as s14, c_value as c14 from hydro_pre_chedh where tag=5.25) t14 on e.id=t14.id and e.sitetype=t14.sitetype
        left join (select id, sitetype, tag, solv as s15, c_value as c15 from hydro_pre_chedh where tag=5.5) t15 on e.id=t15.id and e.sitetype=t15.sitetype
        left join (select id, sitetype, tag, solv as s16, c_value as c16 from hydro_pre_chedh where tag=5.75) t16 on e.id=t16.id and e.sitetype=t16.sitetype
        left join (select id, sitetype, tag, solv as s17, c_value as c17 from hydro_pre_chedh where tag=6) t17 on e.id=t17.id and e.sitetype=t17.sitetype
        left join (select id, sitetype, tag, solv as s18, c_value as c18 from hydro_pre_chedh where tag=6.25) t18 on e.id=t18.id and e.sitetype=t18.sitetype
        left join (select id, sitetype, tag, solv as s19, c_value as c19 from hydro_pre_chedh where tag=6.5) t19 on e.id=t19.id and e.sitetype=t19.sitetype
        left join (select id, sitetype, tag, solv as s20, c_value as c20 from hydro_pre_chedh where tag=6.75) t20 on e.id=t20.id and e.sitetype=t20.sitetype
        left join (select id, sitetype, tag, solv as s21, c_value as c21 from hydro_pre_chedh where tag=7) t21 on e.id=t21.id and e.sitetype=t21.sitetype where c21 is not null;
---- ch
alter table chedh_features add column a_HIS smallint default 0;
alter table chedh_features add column a_CYS smallint default 0;
alter table chedh_features add column b_HIS smallint default 0;
alter table chedh_features add column b_CYS smallint default 0;
alter table chedh_features add column c_HIS smallint default 0;
alter table chedh_features add column c_CYS smallint default 0;
alter table chedh_features add column d_HIS smallint default 0;
alter table chedh_features add column d_CYS smallint default 0;

update chedh_features set a_HIS=1 where resname_a='HIS' and sitetype='ch';
update chedh_features set a_CYS=1 where resname_a='CYS' and sitetype='ch';
update chedh_features set b_HIS=1 where resname_b='HIS' and sitetype='ch';
update chedh_features set b_CYS=1 where resname_b='CYS' and sitetype='ch';
update chedh_features set c_HIS=1 where resname_c='HIS' and sitetype='ch';
update chedh_features set c_CYS=1 where resname_c='CYS' and sitetype='ch';
update chedh_features set d_HIS=1 where resname_d='HIS' and sitetype='ch';
update chedh_features set d_CYS=1 where resname_d='CYS' and sitetype='ch';

alter table chedh_features add column resitype_CC smallint default 0;
alter table chedh_features add column resitype_CH smallint default 0;
alter table chedh_features add column resitype_HH smallint default 0;

update chedh_features set resitype_CC=1 where resi_type='C_C' and sitetype='ch';
update chedh_features set resitype_CH=1 where resi_type='C_H' and sitetype='ch';
update chedh_features set resitype_HH=1 where resi_type='H_H' and sitetype='ch';

---- edh
alter table chedh_features add column a_GLU smallint default 0;
alter table chedh_features add column a_ASP smallint default 0;
alter table chedh_features add column b_GLU smallint default 0;
alter table chedh_features add column b_ASP smallint default 0;
alter table chedh_features add column c_GLU smallint default 0;
alter table chedh_features add column c_ASP smallint default 0;
alter table chedh_features add column d_GLU smallint default 0;
alter table chedh_features add column d_ASP smallint default 0;

update chedh_features set a_HIS=1 where resname_a='HIS' and sitetype='edh';
update chedh_features set a_GLU=1 where resname_a='GLU' and sitetype='edh';
update chedh_features set a_ASP=1 where resname_a='ASP' and sitetype='edh';
update chedh_features set b_HIS=1 where resname_b='HIS' and sitetype='edh';
update chedh_features set b_GLU=1 where resname_b='GLU' and sitetype='edh';
update chedh_features set b_ASP=1 where resname_b='ASP' and sitetype='edh';
update chedh_features set c_HIS=1 where resname_c='HIS' and sitetype='edh';
update chedh_features set c_GLU=1 where resname_c='GLU' and sitetype='edh';
update chedh_features set c_ASP=1 where resname_c='ASP' and sitetype='edh';
update chedh_features set d_HIS=1 where resname_d='HIS' and sitetype='edh';
update chedh_features set d_GLU=1 where resname_d='GLU' and sitetype='edh';
update chedh_features set d_ASP=1 where resname_d='ASP' and sitetype='edh';

alter table chedh_features drop column resname_a;
alter table chedh_features drop column resname_b;
alter table chedh_features drop column resname_c;
alter table chedh_features drop column resname_d;

alter table chedh_features add column resitype_ED smallint default 0;
alter table chedh_features add column resitype_EDH smallint default 0;

update chedh_features set resitype_ED=1 where resi_type='E_D' and sitetype='edh';
update chedh_features set resitype_EDH=1 where resi_type='EDH' and sitetype='edh';
update chedh_features set resitype_HH=1 where resi_type='H_H' and sitetype='edh';

