IRLRealTime
===========

*disclosure, again i did not actually use version control but I am trying to get in better programming practice and starting up.
My real-time codesets which allow you to view fusion code in real-time as well as some corner detection

I still need to put the data up to google drive so that it is backed up a bit better.

--

File		| Description
----------------|--------------
IRLReg.m	| small function which finds the LSR line and returns m, b, and r^2
arr2int.m	| function which turns UDP bytes and flips them into int which can be read by MATLAB
arr2num.m	| generalized version of arr2int.m, flips bytes and returns hex from UDP bytes
calcdisp.m	| calculates the displacement of two 2-d clouds, assuming matches points line up (very shaky)
findclosept.m	| finds the first point in a point cloud close to another.  (DO NOT USE)
findclosepts.m	| finds the closest point in a pointcloud to a given point
findnclusters.m	| finds clusters of points and marks them with their group.  similar to k-means but w/o knowing k.  It's recursive now but in the future I should make it iterative.

--

The only two important files are the ones described below, the other files were me just experimenting with different forms of matching/ feature algorithms.

File			| Description
-------------------------------
udp_cross.m		| Bare minimum for combining LIDAR and vectornav over UDP
udp_line_showoff.m	| Bare minimum for a corner detection over UDP algorithm.


