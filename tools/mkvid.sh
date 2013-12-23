for f in *ppm ; do convert -quality 100 $f `basename $f ppm`jpg; done 
mencoder "mf://*.jpg" -o test.avi -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 