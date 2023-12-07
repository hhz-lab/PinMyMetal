----- for the dist<=2.5:  script  python3 contact_ab.py > contact_ab.csv
drop table contact_ab;
create table contact_ab(
        pdbid char(4),
        cid integer,
        id integer);
\copy contact_ab from './contact_ab.csv' delimiter ','



update zinc_predict_site234_2 t1
        set selected_site=true from
        (select distinct on (pdbid, cid) pdbid, cid, score,id from
                (select a.* , b.score from contact_ab a left join zinc_predict_site234_2 b on a.id=b.id) a
               order by pdbid, cid,score desc) t2
        where t1.id=t2.id;


update zinc_predict_site234_2 t1 set bench='predict' where selected_site is true and bench is null;

