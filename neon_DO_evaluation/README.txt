NEON_DO.csv is generated like so (run on server):

select * from data where upload_id = -900 and variable='DO_mgL' into outfile '/var/lib/mysql-files/NEON_DO.csv' fields terminated by ',' enclosed by '"' lines terminated by '\n';

