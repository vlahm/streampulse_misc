NEON_DO.csv is generated like so (run on server):

select * from data where upload_id = -900 and variable='Discharge_m3s' and DateTime_UTC like '2017%' into outfile '/var/lib/mysql-files/NEON_Q_2017.csv' fields terminated by ',' enclosed by '"' lines terminated by '\n';

