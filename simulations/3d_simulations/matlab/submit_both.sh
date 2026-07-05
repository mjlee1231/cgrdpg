#!/bin/bash
# Submit both n=1000 and n=2000 MATLAB jobs

echo "========================================"
echo "Submitting MATLAB fminunc jobs"
echo "========================================"
echo ""

# Create necessary directories
mkdir -p logs errors results

echo "Submitting n=1000 job..."
job1=$(sbatch submit_matlab_100reps.slurm | awk '{print $4}')
echo "  Job ID: $job1"
echo ""

echo "Submitting n=2000 job..."
job2=$(sbatch submit_matlab_100reps_n2000.slurm | awk '{print $4}')
echo "  Job ID: $job2"
echo ""

echo "========================================"
echo "Both jobs submitted!"
echo "========================================"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "Check logs:"
echo "  tail -f logs/matlab_100reps_${job1}.out      # n=1000"
echo "  tail -f logs/matlab_100reps_n2000_${job2}.out # n=2000"
echo ""
echo "Check results when done:"
echo "  ls -lh results/"
echo ""
