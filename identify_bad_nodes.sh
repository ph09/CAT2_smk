#!/bin/bash
# Script to identify problematic nodes from SLURM log files

echo "Scanning for failed jobs and identifying problematic nodes..."
echo "============================================================"
echo

# Find all error logs with ModuleNotFoundError or ImportError
error_logs=$(find test_install/logs -name "*.err" -type f -exec grep -l "ModuleNotFoundError\|ImportError" {} \;)

if [ -z "$error_logs" ]; then
    echo "No import errors found in logs."
    exit 0
fi

# Extract node names from the corresponding .out files
echo "Failed jobs and their nodes:"
echo "----------------------------"

declare -A bad_nodes

for err_log in $error_logs; do
    # Get corresponding .out file
    out_log="${err_log%.err}.out"
    
    if [ -f "$out_log" ]; then
        # Extract node name from .out file (format: "Node: nodename")
        node=$(grep "^Node:" "$out_log" 2>/dev/null | awk '{print $2}')
        
        if [ -n "$node" ]; then
            echo "  Error log: $err_log"
            echo "  Node: $node"
            echo
            bad_nodes[$node]=1
        fi
    fi
done

if [ ${#bad_nodes[@]} -gt 0 ]; then
    echo "============================================================"
    echo "Problematic nodes identified:"
    echo "----------------------------"
    for node in "${!bad_nodes[@]}"; do
        echo "  - $node"
    done
    
    # Generate exclude string
    exclude_string=$(IFS=,; echo "${!bad_nodes[*]}")
    echo
    echo "============================================================"
    echo "Add this to your input.yaml under default_slurm or specific rules:"
    echo
    echo "  exclude: \"$exclude_string\""
    echo
else
    echo "Could not identify specific nodes from log files."
    echo "Check the .out files manually to see which nodes were used."
fi

