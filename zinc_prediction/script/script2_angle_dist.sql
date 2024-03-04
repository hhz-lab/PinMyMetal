drop table his_zinc_coord23;
create table his_zinc_coord23
        (id smallint,
        zinc_x float,
        zinc_y float,
        zinc_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
	dist_af float,
        dist_bf float,
        dist_cf float);
\copy his_zinc_coord23 from './His_zinc_coord_23.csv' delimiter ','

drop table his_cys_zinc_coord23;
create table his_cys_zinc_coord23
        (id smallint,
        zinc_x float,
        zinc_y float,
        zinc_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
	dist_af float,
        dist_bf float,
        dist_cf float);

\copy his_cys_zinc_coord23 from './cys_his_zinc_coord23.csv' delimiter ','

drop table cys_zinc_coord23;
create table cys_zinc_coord23
        (id smallint,
        zinc_x float,
        zinc_y float,
        zinc_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
	dist_af float,
        dist_bf float,
        dist_cf float
        );
\copy cys_zinc_coord23 from './cys_zinc_coord23.csv' delimiter ','

drop table zinc_coord4;
create table zinc_coord4
        (id smallint,
        zinc_x float,
        zinc_y float,
        zinc_z float,
        atom_a char(4),
        atom_b char(4),
        atom_c char(4),
        atom_d char(4),
	dist_af float,
        dist_bf float,
        dist_cf float,
        dist_df float);
\copy zinc_coord4 from './zinc_coord4.csv' delimiter ','

drop table zinc_predict_site234_2;
create table zinc_predict_site234_2 as
        (select a.*,zinc_x,zinc_y,zinc_z,atom_a,atom_b,atom_c,atom_d,dist_af,dist_bf,dist_cf,dist_df from zinc_predict_site234 a
        left join
        (select * from zinc_coord4 union
        select id,zinc_x,zinc_y,zinc_z,atom_a,atom_b,atom_c,'x' as atom_d,dist_af,dist_bf,dist_cf,-9999 as dist_df from cys_zinc_coord23
        union select id,zinc_x,zinc_y,zinc_z,atom_a,atom_b,atom_c,'x' as atom_d,dist_af,dist_bf,dist_cf,-9999 as dist_df from his_cys_zinc_coord23
        union select id,zinc_x,zinc_y,zinc_z,atom_a,atom_b,atom_c,'x' as atom_d,dist_af,dist_bf,dist_cf,-9999 as dist_df from his_zinc_coord23) b
        on a.id=b.id);

alter table zinc_predict_site234_2 add column select_atom bool default false;
update zinc_predict_site234_2 t1 set select_atom = true
        where atomname_a=atom_a and atomname_b=atom_b and atomname_c=atom_c and atomname_d=atom_d;

delete from zinc_predict_site234_2 where select_atom is false;


--- the dist of CA_a to CA_b

alter table zinc_predict_site234_2 add column CA_dist float;
update zinc_predict_site234_2 set CA_dist = sqrt((CA_x_b-CA_x_a)^2 + (CA_y_b-CA_y_a)^2 + (CA_z_b-CA_z_a)^2);



--- the dist of CB_a to CB_b

alter table zinc_predict_site234_2 add column CB_dist float;
update zinc_predict_site234_2 set CB_dist = sqrt((CB_x_b-CB_x_a)^2 + (CB_y_b-CB_y_a)^2 + (CB_z_b-CB_z_a)^2);



--- the angle of CA to atom
        --- ch 2
alter table zinc_predict_site234_2 add column vector_x_a float;
alter table zinc_predict_site234_2 add column vector_y_a float;
alter table zinc_predict_site234_2 add column vector_z_a float;
update zinc_predict_site234_2 set vector_x_a =x_a-ca_x_a;
update zinc_predict_site234_2 set vector_y_a =y_a-ca_y_a;
update zinc_predict_site234_2 set vector_z_a =z_a-ca_z_a;

alter table zinc_predict_site234_2 add column vector_x_b float;
alter table zinc_predict_site234_2 add column vector_y_b float;
alter table zinc_predict_site234_2 add column vector_z_b float;
update zinc_predict_site234_2 set vector_x_b =x_b-ca_x_b;
update zinc_predict_site234_2 set vector_y_b =y_b-ca_y_b;
update zinc_predict_site234_2 set vector_z_b =z_b-ca_z_b;


alter table zinc_predict_site234_2 add column angle_1 float;
update zinc_predict_site234_2 set angle_1 = (vector_x_a*vector_x_b+vector_y_a*vector_y_b+vector_z_a*vector_z_b);
alter table zinc_predict_site234_2 add column angle_2 float;
update zinc_predict_site234_2 set angle_2 = sqrt(vector_x_a^2 + vector_y_a^2 + vector_z_a^2)*sqrt(vector_x_b^2 + vector_y_b^2 + vector_z_b^2);
alter table zinc_predict_site234_2 add column radian float;
update zinc_predict_site234_2 set radian = ACOS(angle_1/angle_2) where angle_1<> 0 or angle_2 <> 0;
alter table zinc_predict_site234_2 add column angle float;
update zinc_predict_site234_2 set angle = (radian*180)/3.141592653 where angle_1<> 0 or angle_2 <> 0;

alter table zinc_predict_site234_2 drop column vector_x_a;
alter table zinc_predict_site234_2 drop column vector_y_a;
alter table zinc_predict_site234_2 drop column vector_z_a;
alter table zinc_predict_site234_2 drop column vector_x_b;
alter table zinc_predict_site234_2 drop column vector_y_b;
alter table zinc_predict_site234_2 drop column vector_z_b;
alter table zinc_predict_site234_2 drop column angle_1;
alter table zinc_predict_site234_2 drop column angle_2;
alter table zinc_predict_site234_2 drop column radian;

--- the angle1 of CA_a to CB_a and CA_a to CA_b
        --- ch 2
alter table zinc_predict_site234_2 add column vector1_x_a float;
alter table zinc_predict_site234_2 add column vector1_y_a float;
alter table zinc_predict_site234_2 add column vector1_z_a float;
update zinc_predict_site234_2 set vector1_x_a =ca_x_a - cb_x_a;
update zinc_predict_site234_2 set vector1_y_a =ca_y_a - cb_y_a;
update zinc_predict_site234_2 set vector1_z_a =ca_z_a - cb_z_a;

alter table zinc_predict_site234_2 add column vector1_x_b float;
alter table zinc_predict_site234_2 add column vector1_y_b float;
alter table zinc_predict_site234_2 add column vector1_z_b float;
update zinc_predict_site234_2 set vector1_x_b =ca_x_a - ca_x_b;
update zinc_predict_site234_2 set vector1_y_b =ca_y_a - ca_y_b;
update zinc_predict_site234_2 set vector1_z_b =ca_z_a - ca_z_b;


alter table zinc_predict_site234_2 add column angle1_1 float;
update zinc_predict_site234_2 set angle1_1 = (vector1_x_a*vector1_x_b+vector1_y_a*vector1_y_b+vector1_z_a*vector1_z_b);
alter table zinc_predict_site234_2 add column angle1_2 float;
update zinc_predict_site234_2 set angle1_2 = sqrt(vector1_x_a^2 + vector1_y_a^2 + vector1_z_a^2)*sqrt(vector1_x_b^2 + vector1_y_b^2 + vector1_z_b^2);
alter table zinc_predict_site234_2 add column radian1 float;
update zinc_predict_site234_2 set radian1 = ACOS(angle1_1/angle1_2) where angle1_1<> 0 or angle1_2 <> 0;
alter table zinc_predict_site234_2 add column angle1 float;
update zinc_predict_site234_2 set angle1 = (radian1*180)/3.141592653 where angle1_1<> 0 or angle1_2 <> 0;

alter table zinc_predict_site234_2 drop column vector1_x_a;
alter table zinc_predict_site234_2 drop column vector1_y_a;
alter table zinc_predict_site234_2 drop column vector1_z_a;
alter table zinc_predict_site234_2 drop column vector1_x_b;
alter table zinc_predict_site234_2 drop column vector1_y_b;
alter table zinc_predict_site234_2 drop column vector1_z_b;
alter table zinc_predict_site234_2 drop column angle1_1;
alter table zinc_predict_site234_2 drop column angle1_2;
alter table zinc_predict_site234_2 drop column radian1;

--- the angle2 of CA_B to CB_B and CA_B to CA_a
alter table zinc_predict_site234_2 add column vector2_x_a float;
alter table zinc_predict_site234_2 add column vector2_y_a float;
alter table zinc_predict_site234_2 add column vector2_z_a float;
update zinc_predict_site234_2 set vector2_x_a =ca_x_b - cb_x_b;
update zinc_predict_site234_2 set vector2_y_a =ca_y_b - cb_y_b;
update zinc_predict_site234_2 set vector2_z_a =ca_z_b - cb_z_b;

alter table zinc_predict_site234_2 add column vector2_x_b float;
alter table zinc_predict_site234_2 add column vector2_y_b float;
alter table zinc_predict_site234_2 add column vector2_z_b float;
update zinc_predict_site234_2 set vector2_x_b =ca_x_b - ca_x_a;
update zinc_predict_site234_2 set vector2_y_b =ca_y_b - ca_y_a;
update zinc_predict_site234_2 set vector2_z_b =ca_z_b - ca_z_a;


alter table zinc_predict_site234_2 add column angle2_1 float;
update zinc_predict_site234_2 set angle2_1 = (vector2_x_a*vector2_x_b+vector2_y_a*vector2_y_b+vector2_z_a*vector2_z_b);
alter table zinc_predict_site234_2 add column angle2_2 float;
update zinc_predict_site234_2 set angle2_2 = sqrt(vector2_x_a^2 + vector2_y_a^2 + vector2_z_a^2)*sqrt(vector2_x_b^2 + vector2_y_b^2 + vector2_z_b^2);
alter table zinc_predict_site234_2 add column radian2 float;
update zinc_predict_site234_2 set radian2 = ACOS(angle2_1/angle2_2) where angle2_1<> 0 or angle2_2 <> 0;
alter table zinc_predict_site234_2 add column angle2 float;
update zinc_predict_site234_2 set angle2 = (radian2*180)/3.141592653 where angle2_1<> 0 or angle2_2 <> 0;

alter table zinc_predict_site234_2 drop column vector2_x_a;
alter table zinc_predict_site234_2 drop column vector2_y_a;
alter table zinc_predict_site234_2 drop column vector2_z_a;
alter table zinc_predict_site234_2 drop column vector2_x_b;
alter table zinc_predict_site234_2 drop column vector2_y_b;
alter table zinc_predict_site234_2 drop column vector2_z_b;
alter table zinc_predict_site234_2 drop column angle2_1;
alter table zinc_predict_site234_2 drop column angle2_2;
alter table zinc_predict_site234_2 drop column radian2;

--- the angle3 of CA_a to CB_a and CA_b to CB_b

        --- ch 2
alter table zinc_predict_site234_2 add column vector3_x_a float;
alter table zinc_predict_site234_2 add column vector3_y_a float;
alter table zinc_predict_site234_2 add column vector3_z_a float;
update zinc_predict_site234_2 set vector3_x_a =ca_x_a - cb_x_a;
update zinc_predict_site234_2 set vector3_y_a =ca_y_a - cb_y_a;
update zinc_predict_site234_2 set vector3_z_a =ca_z_a - cb_z_a;

alter table zinc_predict_site234_2 add column vector3_x_b float;
alter table zinc_predict_site234_2 add column vector3_y_b float;
alter table zinc_predict_site234_2 add column vector3_z_b float;
update zinc_predict_site234_2 set vector3_x_b =ca_x_b - cb_x_b;
update zinc_predict_site234_2 set vector3_y_b =ca_y_b - cb_y_b;
update zinc_predict_site234_2 set vector3_z_b =ca_z_b - cb_z_b;


alter table zinc_predict_site234_2 add column angle3_1 float;
update zinc_predict_site234_2 set angle3_1 = (vector3_x_a*vector3_x_b+vector3_y_a*vector3_y_b+vector3_z_a*vector3_z_b);
alter table zinc_predict_site234_2 add column angle3_2 float;
update zinc_predict_site234_2 set angle3_2 = sqrt(vector3_x_a^2 + vector3_y_a^2 + vector3_z_a^2)*sqrt(vector3_x_b^2 + vector3_y_b^2 + vector3_z_b^2);
alter table zinc_predict_site234_2 add column radian3 float;
update zinc_predict_site234_2 set radian3 = ACOS(angle3_1/angle3_2) where angle3_1<> 0 or angle3_2 <> 0;
alter table zinc_predict_site234_2 add column angle3 float;
update zinc_predict_site234_2 set angle3 = (radian3*180)/3.141592653 where angle3_1<> 0 or angle3_2 <> 0;

alter table zinc_predict_site234_2 drop column vector3_x_a ;
alter table zinc_predict_site234_2 drop column vector3_y_a ;
alter table zinc_predict_site234_2 drop column vector3_z_a ;
alter table zinc_predict_site234_2 drop column vector3_x_b ;
alter table zinc_predict_site234_2 drop column vector3_y_b ;
alter table zinc_predict_site234_2 drop column vector3_z_b ;
alter table zinc_predict_site234_2 drop column angle3_1 ;
alter table zinc_predict_site234_2 drop column angle3_2 ;
alter table zinc_predict_site234_2 drop column radian3 ;
--- the score of dist and angle


alter table zinc_predict_site234_2 add column atom_type text;
update zinc_predict_site234_2 set atom_type = 'S_N' where (atomname_a = '-SG-' and atomname_b = '-NE2' )
        or (atomname_b = '-SG-' and atomname_a = '-NE2')
        or (atomname_a = '-SG-' and atomname_b = '-ND1' )
        or (atomname_b = '-SG-' and atomname_a = '-ND1');

update zinc_predict_site234_2 set atom_type = 'N_N' where (atomname_a = '-NE2' and atomname_b = '-ND1' )
        or (atomname_b = '-NE2' and atomname_a = '-ND1')
        or (atomname_a = '-ND1' and atomname_b = '-ND1')
        or (atomname_a = '-NE2' and atomname_b = '-NE2' );

update zinc_predict_site234_2 set atom_type = 'S_S' where (atomname_a = '-SG-' and atomname_b = '-SG-' );

alter table zinc_predict_site234_2 add column atom_dist_2 float;
alter table zinc_predict_site234_2 add column ca_dist_2 float;
alter table zinc_predict_site234_2 add column cb_dist_2 float;
alter table zinc_predict_site234_2 add column angle_2 float;
alter table zinc_predict_site234_2 add column angle1_2 float;
alter table zinc_predict_site234_2 add column angle2_2 float;
alter table zinc_predict_site234_2 add column angle3_2 float;

update zinc_predict_site234_2 set atom_dist_2 = cast(dist_ab/0.5 as integer)*0.5 ;
update zinc_predict_site234_2 set ca_dist_2 = cast(CA_dist/2 as integer)*2 ;
update zinc_predict_site234_2 set cb_dist_2 = cast(CB_dist/2 as integer)*2 ;
update zinc_predict_site234_2 set angle_2 = cast(angle/20 as integer)*20;
update zinc_predict_site234_2 set angle1_2 = cast(angle1/20 as integer)*20;
update zinc_predict_site234_2 set angle2_2 = cast(angle2/20 as integer)*20;
update zinc_predict_site234_2 set angle3_2 = cast(angle3/20 as integer)*20;


alter table zinc_predict_site234_2 add column atom_dist_score numeric default 0;
alter table zinc_predict_site234_2 add column ca_dist_score numeric default 0;
alter table zinc_predict_site234_2 add column cb_dist_score numeric default 0;
alter table zinc_predict_site234_2 add column angle_score numeric default 0;
alter table zinc_predict_site234_2 add column angle1_score numeric default 0;
alter table zinc_predict_site234_2 add column angle2_score numeric default 0;
alter table zinc_predict_site234_2 add column angle3_score numeric default 0;
alter table zinc_predict_site234_2 add column  score numeric default 0;

update zinc_predict_site234_2 t1
        set atom_dist_score = t2.value_score from (select * from dist_angle_score ) t2
        where t1.atom_dist_2 = t2.value and t1.atom_type=t2.atom_type and value_type = 'atom_dist' ;
update zinc_predict_site234_2 t1
        set ca_dist_score = t2.value_score from (select * from dist_angle_score ) t2
        where t1.ca_dist_2 = t2.value and t1.atom_type=t2.atom_type and value_type = 'ca_dist' ;
update zinc_predict_site234_2 t1
        set cb_dist_score = t2.value_score from (select * from dist_angle_score ) t2
        where t1.cb_dist_2 = t2.value and t1.atom_type=t2.atom_type and value_type = 'cb_dist' ;
update zinc_predict_site234_2 t1
        set angle_score = t2.value_score from (select * from dist_angle_score ) t2
        where t1.angle_2 = t2.value and t1.atom_type=t2.atom_type and value_type = 'angle' ;
update zinc_predict_site234_2 t1
        set angle1_score = t2.value_score from (select * from dist_angle_score ) t2
        where t1.angle1_2 = t2.value and t1.atom_type=t2.atom_type and value_type = 'angle1' ;
update zinc_predict_site234_2 t1
        set angle2_score = t2.value_score from (select * from dist_angle_score ) t2
        where t1.angle2_2 = t2.value and t1.atom_type=t2.atom_type and value_type = 'angle2' ;
update zinc_predict_site234_2 t1
        set angle3_score = t2.value_score from (select * from dist_angle_score ) t2
        where t1.angle3_2 = t2.value and t1.atom_type=t2.atom_type and value_type = 'angle3' ;

update zinc_predict_site234_2 set
        score = cast((atom_dist_score+ca_dist_score+cb_dist_score+angle_score+angle1_score+angle2_score+angle3_score) as decimal(18,3));

