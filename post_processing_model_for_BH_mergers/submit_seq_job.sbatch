#!/bin/bash/

desired_values=(1 5 25 125 625 3125)
# Iterate through the values in the array
for i in "${desired_values[@]}"; do
   echo "Welcome $i times"
   sbatch --dependency=singleton submit_job.sh $i
done
