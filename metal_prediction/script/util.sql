-- FUNCTION AGGREGATING TEXT FIELDS WITH COMMAS
create language plpgsql;
drop function comma_aggregate(text,text) cascade;
create function comma_aggregate(text,text) returns text as '
begin
  if (length($2) > 0 ) then
    if (length($1) >0 ) then
      return $1 || '', '' || $2;
    else
      return $2;
    end if;
  else
    return $1;
  end if;
end;
' language 'plpgsql';


-- create the aggregate function
create aggregate conc_comma (
	basetype=text,
	sfunc=comma_aggregate,
        stype=text,
	initcond=''
);


-- FUNCTION AGGREGATING TEXT FIELDS WITH COMMAS
create language plpgsql;
drop function semicolon_aggregate(text,text) cascade;
create function semicolon_aggregate(text,text) returns text as '
begin
  if (length($2) > 0 ) then
    if (length($1) >0 ) then
      return $1 || ''; '' || $2;
    else
      return $2;
    end if;
  else
    return $1;
  end if;
end;
' language 'plpgsql';


-- create the aggregate function
create aggregate conc_semicolon (
	basetype=text,
	sfunc=semicolon_aggregate,
        stype=text,
	initcond=''
);


-- FUNCTION AGGREGATING NUMERIC FIELDS TO GIVE MEDIAN VALUE
--CREATE AGGREGATE median(numeric) (
--  SFUNC=array_append,
--  STYPE=numeric[],
--  FINALFUNC=array_median
--);
--
--CREATE AGGREGATE median(timestamp) (
--  SFUNC=array_append,
--  STYPE=timestamp[],
--  FINALFUNC=array_median
--);

CREATE OR REPLACE FUNCTION array_median(numeric[])
  RETURNS numeric AS
$$
    SELECT CASE WHEN array_upper($1,1) = 0 THEN null ELSE asorted[ceiling(array_upper(asorted,1)/2.0)] END
    FROM (SELECT ARRAY(SELECT ($1)[n] FROM
generate_series(1, array_upper($1, 1)) AS n
    WHERE ($1)[n] IS NOT NULL
            ORDER BY ($1)[n]
) As asorted) As foo ;
$$
  LANGUAGE 'sql' IMMUTABLE;

CREATE OR REPLACE FUNCTION array_median(timestamp[])
  RETURNS timestamp AS
$$
    SELECT CASE WHEN array_upper($1,1) = 0 THEN null ELSE asorted[ceiling(array_upper(asorted,1)/2.0)] END
    FROM (SELECT ARRAY(SELECT ($1)[n] FROM
generate_series(1, array_upper($1, 1)) AS n
    WHERE ($1)[n] IS NOT NULL
            ORDER BY ($1)[n]
) As asorted) As foo ;
$$
  LANGUAGE 'sql' IMMUTABLE;


-- create the aggregate function
CREATE AGGREGATE median(numeric) (
  SFUNC=array_append,
  STYPE=numeric[],
  FINALFUNC=array_median
);

CREATE AGGREGATE median(timestamp) (
  SFUNC=array_append,
  STYPE=timestamp[],
  FINALFUNC=array_median
);

