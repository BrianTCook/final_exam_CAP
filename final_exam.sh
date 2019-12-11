#'should be run in ~/amuse directory'
#'should put into ~/final_exam_CAP for safekeeping'

echo 'run script'
./amuse.sh ~/Desktop/comp_astro/final_exam_CAP/source.py

echo 'convert collection of .png files into .gif'
convert -delay 10 'frame_*.png' -loop 0 cluster_evolution.gif
rm -rf frame_*.png

mv *.png ~/Desktop/comp_astro/final_exam_CAP/
mv *.gif ~/Desktop/comp_astro/final_exam_CAP/

echo 'copy shell script into same directory for safekeeping'
cp final_exam.sh ~/Desktop/comp_astro/final_exam_CAP/
