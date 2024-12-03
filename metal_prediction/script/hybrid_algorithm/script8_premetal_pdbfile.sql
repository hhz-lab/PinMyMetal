--- pre metals info
drop table presites_info_neighbor;
create table presites_info_neighbor
        (id real,
	sitetype char(3),
        pdbid char(4),
        chainid_ion char(4),
        resseq_ion integer,
        residueid_ion integer);
\copy presites_info_neighbor from './pre_metals_info.csv' delimiter ','

