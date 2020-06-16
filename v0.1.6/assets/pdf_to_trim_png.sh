muhfile="full_feedback_control_system"
convert -density 300 $muhfile.pdf $muhfile.png
convert $muhfile.png -trim $muhfile.png
