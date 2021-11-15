# 2-D boundary description
# "#" is the comment line signal
# *SIZE= How many loops do we have in the boundary :
# go ccw for externl and CW for internal boundaries 
#*LOOP is to list the segents in a loops 
# Dont put comments in the loop section 
#each segement has 
#Number of points, Order of spline, 
#whther it is periodic or not (for smooth closed loops like holes) 
#points x,y,x,y,x,y etc 
#*END is the end of the file 
*SIZE
2
*LOOP
20,1, 0,-1,-1, 1,-1
20,1,0, 1,-1, 1,1
20,1,0, 1,1, -1,1
20,1,0, -1,1, -1,-1
*LOOP
40,3, 1,  0,-0.5, -0.5,0, 0,0.5, 0.5,0, 0,-0.5
*END