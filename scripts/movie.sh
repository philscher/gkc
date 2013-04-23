mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc copy -oac copy -o output.avi
mencoder mf://plot*.png -fps 5 -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=2400 -ffourcc DIVX -vf scale=640:480 -o LinearETG.avi
