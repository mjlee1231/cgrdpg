#!/bin/bash
# Quick test of surrogate implementation on HPC (no SLURM submission)
# Run directly: bash test_surrogate_hpc.sh

echo "========================================"
echo "Testing Surrogate Likelihood on HPC"
echo "========================================"
echo ""

module load matlab

echo "Running quick test (n=500)..."
matlab -nodisplay -nosplash -nodesktop -r "quick_test_surrogate; exit;"

exit_code=$?

echo ""
echo "========================================"
if [ $exit_code -eq 0 ]; then
    echo "Test completed successfully!"
    echo "Check output above for improvement %"
    echo ""
    echo "If improvement >50%, the surrogate works!"
    echo "Then we can update the 100-rep scripts."
else
    echo "Test failed with exit code: $exit_code"
    echo "Check for errors above."
fi
echo "========================================"
