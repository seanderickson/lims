/*** TODO: this script is provisional ** may not be necessary ** 
*** keep as an example **
***/

/** ecommons id's must be null or unique **/
/** identify duplicates: select count(*), ecommons_id, array_agg(screensaver_user_id) from screensaver_user where ecommons_id is not null group by ecommons_id order by count desc; **/

/** ecommons id's either must be null or unique.  This is an old user account for maria chmurra **/
update screensaver_user set ecommons_id =null where screensaver_user_id = 129;

/**          2 | jga10       | {1408,3166} 3166 |     244 | 2008-02-13 00:00:00 | John       | Gerald , 1408 |     325 | 2007-10-29 00:00:00 | John       | Albeck **/
update screensaver_user set ecommons_id =null where screensaver_user_id = 1408;

/**          2 | me44        | {1234,866}   866 |      14 | 2007-05-24 00:00:00 | Ernebjerg  | Morten , 1234 |     331 | 2007-05-24 00:00:00 | Morten     | Ernebjerg **/
update screensaver_user set ecommons_id =null where screensaver_user_id = 866;

/**          2 | rcs12       | {4027,830}  830 |     346 | 2007-04-25 00:00:00     | Ruchir     | Shah **/
update screensaver_user set ecommons_id =null where screensaver_user_id = 830;
update screensaver_user set ecommons_id =null where screensaver_user_id = 4027;


