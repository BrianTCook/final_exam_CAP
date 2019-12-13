#'should be run in ~/amuse directory'
#'should put into ~/cap_assignments/final_project for safekeeping'

echo 'run script'
./amuse.sh ~/Desktop/comp_astro/final_exam_CAP/run_Neptune_with_disk_and_black_hole.py

cp *.png ~/Desktop/comp_astro/final_exam_CAP/

echo 'copy shell script into same directory for safekeeping'
cp final_exam.sh ~/Desktop/comp_astro/final_exam_CAP/
