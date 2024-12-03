-------chem features
drop table chemfeatures_pre_ch;
create table chemfeatures_pre_ch
        (id integer, bin1_class0 real, bin1_class1 real, bin1_class2 real, bin1_class3 real, bin1_class4 real, bin1_class5 real, bin1_class6 real, bin1_class7 real, bin2_class0 real, bin2_class1 real, bin2_class2 real, bin2_class3 real, bin2_class4 real, bin2_class5 real, bin2_class6 real, bin2_class7 real, bin3_class0 real, bin3_class1 real, bin3_class2 real, bin3_class3 real, bin3_class4 real, bin3_class5 real, bin3_class6 real, bin3_class7 real);
\copy chemfeatures_pre_ch from './chem_ch.csv' delimiter ',' csv header

drop table chemfeatures_pre_edh;
create table chemfeatures_pre_edh
        (id integer, bin1_class0 real, bin1_class1 real, bin1_class2 real, bin1_class3 real, bin1_class4 real, bin1_class5 real, bin1_class6 real, bin1_class7 real, bin2_class0 real, bin2_class1 real, bin2_class2 real, bin2_class3 real, bin2_class4 real, bin2_class5 real, bin2_class6 real, bin2_class7 real, bin3_class0 real, bin3_class1 real, bin3_class2 real, bin3_class3 real, bin3_class4 real, bin3_class5 real, bin3_class6 real, bin3_class7 real);
\copy chemfeatures_pre_edh from './chem_edh.csv' delimiter ',' csv header

------ ionbinding features
drop table ion_bindingsites;
CREATE TABLE ion_bindingsites (
    pdbfileid           integer,
    bindingsiteid       integer NOT NULL,
    chainid             char(3),
    residueid_ion       smallint NOT NULL,
    resname_ion         char(3),
    atomid_ion          integer NOT NULL,
    atomname_ion        char(4),
    protons_ion         smallint,
    oxidation_state     smallint,
    bfactor_ion         real,
    bfactor_env_avg     real,
    occupancy_ion       real,
    occupancy_env_avg   real,
    coordnum_inner      smallint,
    coordnum_bidentate_merged smallint,
    geometry_type       smallint,
    geometry_bidentate  smallint,       -- number of bidentate-pseudo-atom used
    geometry_pseudo     smallint,       -- number of vertex missed
    geometry_distort    boolean,        -- RMSD angle bigger than 15
    rmsd_geom_angle     real,           -- RMSD angle
    num_bidentate_all   smallint,
    distance_avg        real,
    distance_min        real,
    distance_max        real,
    valence_3a          real,
    vecsum_3a           real,
    coord_number_4a     smallint,
    num_metal_4a        smallint,
    valence_4a          real,
    vecsum_4a           real,
    bvs_sodium          real,
    bvs_magnesium       real,
    bvs_potassium       real,
    bvs_calcium         real,
    bvs_manganese       real,
    bvs_iron            real,
    bvs_cobalt          real,
    bvs_nickel          real,
    bvs_copper          real,
    bvs_zinc            real
);
\copy ion_bindingsites from './presite_feature_ion.txt' delimiter '|'


drop table ionresidueid_id;
create table ionresidueid_id (
        id real,
        sitetype char(3),
        pdbid char(4),
        chainid_ion char(3),
        resseq_ion integer,
        pdbfileid integer,
        residueid_ion smallint);
\copy ionresidueid_id from './presite_feature_info.txt' delimiter '|'

-------6. merge features
  --------- CH
drop table classmodel_feature_neighbor_ch;
create table classmodel_feature_neighbor_ch as
        select t1.*,
         bin1_class0,bin1_class1,bin1_class2,bin1_class3,bin1_class4,bin1_class5,bin1_class6,bin1_class7,
        bin2_class0,bin2_class1,bin2_class2,bin2_class3,bin2_class4,bin2_class5,bin2_class6,bin2_class7,
        bin3_class0,bin3_class1,bin3_class2,bin3_class3,bin3_class4,bin3_class5,bin3_class6,bin3_class7,
        coordnum_inner,coordnum_bidentate_merged,geometry_type,geometry_bidentate,geometry_pseudo,geometry_distort,rmsd_geom_angle,num_bidentate_all,distance_avg,distance_min,distance_max,valence_3a,vecsum_3a,coord_number_4a,num_metal_4a,valence_4a,vecsum_4a,bvs_sodium,bvs_magnesium,bvs_potassium,bvs_calcium,bvs_manganese,bvs_iron,bvs_cobalt,bvs_nickel,bvs_copper,bvs_zinc
        from (
        select a.pdbfileid,a.residueid_ion,b.* from (select * from ionresidueid_id where sitetype='ch') a left join
        (select id, sitetype, site_count,
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21
	from chedh_features where sitetype='ch') b
        on a.id = b.id ) t1
        left join ion_bindingsites t2
        on t1.pdbfileid=t2.pdbfileid and t1.residueid_ion=t2.residueid_ion
        left join chemfeatures_pre_ch t3
        on t1.id=t3.id;

alter table classmodel_feature_neighbor_ch add column geometry_distort2 smallint;
update classmodel_feature_neighbor_ch set geometry_distort2 = 1 where geometry_distort is true;
update classmodel_feature_neighbor_ch set geometry_distort2 = 0 where geometry_distort is false;

alter table classmodel_feature_neighbor_ch drop column geometry_distort;
alter table classmodel_feature_neighbor_ch rename geometry_distort2 to geometry_distort;

         ----------- EDH

drop table classmodel_feature_neighbor_edh;
create table classmodel_feature_neighbor_edh as
        select t1.*,
         bin1_class0,bin1_class1,bin1_class2,bin1_class3,bin1_class4,bin1_class5,bin1_class6,bin1_class7,
        bin2_class0,bin2_class1,bin2_class2,bin2_class3,bin2_class4,bin2_class5,bin2_class6,bin2_class7,
        bin3_class0,bin3_class1,bin3_class2,bin3_class3,bin3_class4,bin3_class5,bin3_class6,bin3_class7,
        coordnum_inner,coordnum_bidentate_merged,geometry_type,geometry_bidentate,geometry_pseudo,geometry_distort,rmsd_geom_angle,num_bidentate_all,distance_avg,distance_min,distance_max,valence_3a,vecsum_3a,coord_number_4a,num_metal_4a,valence_4a,vecsum_4a,bvs_sodium,bvs_magnesium,bvs_potassium,bvs_calcium,bvs_manganese,bvs_iron,bvs_cobalt,bvs_nickel,bvs_copper,bvs_zinc
        from (
        select a.pdbfileid,a.residueid_ion,b.* from (select * from ionresidueid_id where sitetype='edh') a left join
        (select id,sitetype, site_count,
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,
    a_his,a_glu,a_asp,b_his,b_glu,b_asp,c_his,c_glu,c_asp,d_his,d_glu,d_asp,resitype_ed,resitype_edh,resitype_hh 
	from chedh_features where sitetype='edh') b
        on a.id = b.id ) t1
        left join ion_bindingsites t2
        on t1.pdbfileid=t2.pdbfileid and t1.residueid_ion=t2.residueid_ion
        left join chemfeatures_pre_edh t3
        on t1.id=t3.id;

alter table classmodel_feature_neighbor_edh add column geometry_distort2 smallint;
update classmodel_feature_neighbor_edh set geometry_distort2 = 1 where geometry_distort is true;
update classmodel_feature_neighbor_edh set geometry_distort2 = 0 where geometry_distort is false;

alter table classmodel_feature_neighbor_edh drop column geometry_distort;
alter table classmodel_feature_neighbor_edh rename geometry_distort2 to geometry_distort;

--- merge ch and edh sites
drop table classmodel_feature_neighbor_chedh;
create table classmodel_feature_neighbor_chedh as
        select distinct id, sitetype,
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,
        bin1_class0,bin1_class1,bin1_class2,bin1_class3,bin1_class4,bin1_class5,bin1_class6,bin1_class7,
        bin2_class0,bin2_class1,bin2_class2,bin2_class3,bin2_class4,bin2_class5,bin2_class6,bin2_class7,
        bin3_class0,bin3_class1,bin3_class2,bin3_class3,bin3_class4,bin3_class5,bin3_class6,bin3_class7,
        coordnum_inner,coordnum_bidentate_merged,geometry_type,geometry_bidentate,geometry_pseudo,geometry_distort,rmsd_geom_angle,num_bidentate_all,distance_avg,distance_min,distance_max,valence_3a,vecsum_3a,coord_number_4a,num_metal_4a,valence_4a,vecsum_4a,bvs_sodium,bvs_magnesium,bvs_potassium,bvs_calcium,bvs_manganese,bvs_iron,bvs_cobalt,bvs_nickel,bvs_copper,bvs_zinc from classmodel_feature_neighbor_ch union
        select distinct id, sitetype,
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
        s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,
        bin1_class0,bin1_class1,bin1_class2,bin1_class3,bin1_class4,bin1_class5,bin1_class6,bin1_class7,
        bin2_class0,bin2_class1,bin2_class2,bin2_class3,bin2_class4,bin2_class5,bin2_class6,bin2_class7,
        bin3_class0,bin3_class1,bin3_class2,bin3_class3,bin3_class4,bin3_class5,bin3_class6,bin3_class7,
        coordnum_inner,coordnum_bidentate_merged,geometry_type,geometry_bidentate,geometry_pseudo,geometry_distort,rmsd_geom_angle,num_bidentate_all,distance_avg,distance_min,distance_max,valence_3a,vecsum_3a,coord_number_4a,num_metal_4a,valence_4a,vecsum_4a,bvs_sodium,bvs_magnesium,bvs_potassium,bvs_calcium,bvs_manganese,bvs_iron,bvs_cobalt,bvs_nickel,bvs_copper,bvs_zinc from classmodel_feature_neighbor_edh;


