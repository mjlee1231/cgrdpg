#!/bin/bash
# Check if MATLAB is available on HPC and show available versions

echo "========================================"
echo "Checking MATLAB availability on HPC"
echo "========================================"
echo ""

# Check if module command exists
if command -v module &> /dev/null; then
    echo "✓ Module system available"
    echo ""

    # Search for MATLAB modules
    echo "Available MATLAB modules:"
    module avail matlab 2>&1 | grep -i matlab || echo "  No MATLAB modules found"
    echo ""

    # Try to load MATLAB
    echo "Attempting to load MATLAB module..."
    module load matlab 2>&1

    if command -v matlab &> /dev/null; then
        echo "✓ MATLAB loaded successfully!"
        echo ""
        matlab -batch "version; ver('parallel')"
    else
        echo "✗ MATLAB not found after module load"
        echo ""
        echo "Try one of these:"
        echo "  module load matlab/R2023a"
        echo "  module load matlab/R2022b"
        echo "  module spider matlab"
    fi
else
    echo "✗ Module system not available"
    echo ""

    # Check if MATLAB is in PATH
    if command -v matlab &> /dev/null; then
        echo "✓ MATLAB found in PATH"
        matlab -batch "version"
    else
        echo "✗ MATLAB not found"
    fi
fi

echo ""
echo "========================================"
echo "System information:"
echo "========================================"
echo "Hostname: $(hostname)"
echo "CPUs: $(nproc)"
echo "Memory: $(free -h | grep Mem | awk '{print $2}')"
echo ""
