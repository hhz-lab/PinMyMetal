-------- CH
drop table ch_metal_coord23;
create table ch_metal_coord23
        (id integer,
        metal_x float,
        metal_y float,
        metal_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
	dist_am float,
	dist_bm float,
	dist_cm float);

\copy ch_metal_coord23 from './CH_metal_coord23.csv'

drop table ch_metal_coord4;
create table ch_metal_coord4
        (id integer,
        metal_x float,
        metal_y float,
        metal_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
        atom_d char(4),
	dist_am float,
        dist_bm float,
        dist_cm float,
	dist_dm float);
\copy ch_metal_coord4 from './CH_metal_coord4.csv' delimiter ','


drop table zncu_predict_sites_2;
create table zncu_predict_sites_2 as
        select a.*,metal_x,metal_y,metal_z,atom_a,atom_b,atom_c,atom_d,dist_am,dist_bm,dist_cm,dist_dm from
                (select id,metal_x,metal_y,metal_z,atom_a,atom_b,atom_c,atom_d,dist_am,dist_bm,dist_cm,dist_dm from ch_metal_coord4
        union select id,metal_x,metal_y,metal_z,atom_a,atom_b,atom_c,'NaN' as atom_d,dist_am,dist_bm,dist_cm, 0 as dist_dm from ch_metal_coord23) b
        left join
                (select distinct id,pdbid,pdbfileid,residueid_a,residueid_b,residueid_c,residueid_d,
                        resname_a,resname_b,resname_c,resname_d,
                        resseq_a,resseq_b,resseq_c,resseq_d,
                        chainid_a,chainid_b,chainid_c,chainid_d,
                        site_count,resi_type,
                        contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd from
                        zncu_predict_sites) a
        on a.id=b.id;




--------- EDH
drop table edh_metal_coord;
create table edh_metal_coord
        (id integer,
        metal_x float,
        metal_y float,
        metal_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
        atom_d char(4),
        atom_e char(4),
	dist_am float,
        dist_bm float,
        dist_cm float,
        dist_dm float,
	dist_em float);

\copy edh_metal_coord from './EDH_metal_coord2345.csv'

drop table metal_coord6;
create table metal_coord6
        (id integer,
        metal_x float,
        metal_y float,
        metal_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
        atom_d char(4),
        atom_e char(4),
        atom_f char(4),
	dist_am float,
        dist_bm float,
        dist_cm float,
        dist_dm float,
        dist_em float,
	dist_fm float);
\copy metal_coord6 from './EDH_metal_coord6.csv' delimiter ','

drop table metal_predict_sites_2;
create table metal_predict_sites_2 as
        select a.*,metal_x,metal_y,metal_z,atom_a,atom_b,atom_c,atom_d,atom_e,atom_f,dist_am,dist_bm,dist_cm,dist_dm,dist_em,dist_fm
        from (select id,metal_x,metal_y,metal_z,atom_a,atom_b,atom_c,atom_d,atom_e,atom_f,dist_am,dist_bm,dist_cm,dist_dm,dist_em,dist_fm from metal_coord6
        union select id,metal_x,metal_y,metal_z,atom_a,atom_b,atom_c,atom_d,atom_e,'NaN' as atom_f,dist_am,dist_bm,dist_cm,dist_dm,dist_em,0 as dist_fm from edh_metal_coord) b
        left join
                (select distinct id,pdbid,pdbfileid,residueid_a,residueid_b,residueid_c,residueid_d,residueid_e,residueid_f,
                        resname_a,resname_b,resname_c,resname_d,resname_e,resname_f,
                        resseq_a,resseq_b,resseq_c,resseq_d,resseq_e,resseq_f,
                        chainid_a,chainid_b,chainid_c,chainid_d,chainid_e,chainid_f,
                        site_count,resi_type,
                        contact_ab,contact_bc,contact_ca,contact_ad,contact_bd,contact_cd,contact_ae,contact_be,contact_ce,contact_de,contact_af,contact_bf,contact_cf,contact_df,contact_ef from
                        metal_predict_sites) a
        on a.id=b.id;



