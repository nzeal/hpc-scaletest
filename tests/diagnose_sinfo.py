#!/usr/bin/env python3
"""
Diagnostic script to debug sinfo partition detection.
Run this on the HPC system to see what sinfo actually returns.
"""

import subprocess
import sys

def run_sinfo_variants(partition):
    """Try different sinfo command variants to see which works."""
    
    commands = [
        # Current framework command
        ['sinfo', '-p', partition, '-o', '%D %c %m %G', '-h'],
        
        # Alternative formats
        ['sinfo', '-p', partition, '-o', '%P %D %c %m %G', '-h'],
        ['sinfo', '-p', partition, '-o', '%all', '-h'],
        ['sinfo', '-p', partition, '--Format=partition,nodes,cpus,memory,gres', '-h'],
        ['sinfo', '-p', partition, '--Format=gres:30', '-h'],
        
        # Without -h (with header)
        ['sinfo', '-p', partition, '-o', '%D %c %m %G'],
        
        # Node-level query
        ['sinfo', '-p', partition, '-N', '-o', '%N %c %m %G', '-h'],
    ]
    
    print(f"\n{'='*80}")
    print(f"SINFO DIAGNOSTIC FOR PARTITION: {partition}")
    print(f"{'='*80}\n")
    
    for i, cmd in enumerate(commands, 1):
        print(f"\n[{i}] Command: {' '.join(cmd)}")
        print("-" * 80)
        
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            
            print(f"Return code: {result.returncode}")
            
            if result.returncode == 0:
                if result.stdout.strip():
                    print(f"STDOUT (length={len(result.stdout)}):")
                    print(result.stdout)
                else:
                    print("STDOUT: <empty>")
            else:
                print(f"FAILED with return code {result.returncode}")
            
            if result.stderr.strip():
                print(f"STDERR:")
                print(result.stderr)
                
        except subprocess.TimeoutExpired:
            print("TIMEOUT after 10 seconds")
        except Exception as e:
            print(f"ERROR: {e}")
    
    # Also try querying node directly
    print(f"\n{'='*80}")
    print("TRYING NODE-LEVEL QUERY")
    print(f"{'='*80}\n")
    
    try:
        # Get a node from this partition
        node_cmd = ['sinfo', '-p', partition, '-N', '-h', '-o', '%N']
        result = subprocess.run(node_cmd, capture_output=True, text=True, timeout=10)
        
        if result.returncode == 0 and result.stdout.strip():
            nodes = result.stdout.strip().split('\n')
            first_node = nodes[0].strip()
            
            print(f"First node in partition: {first_node}")
            print("-" * 80)
            
            # Query this specific node
            node_query = ['sinfo', '-n', first_node, '-o', '%N %c %m %G', '-h']
            print(f"Command: {' '.join(node_query)}")
            
            result = subprocess.run(node_query, capture_output=True, text=True, timeout=10)
            print(f"Return code: {result.returncode}")
            print(f"Output:\n{result.stdout}")
            
    except Exception as e:
        print(f"Node query failed: {e}")

    # Try scontrol show partition
    print(f"\n{'='*80}")
    print("SCONTROL SHOW PARTITION")
    print(f"{'='*80}\n")
    
    try:
        scontrol_cmd = ['scontrol', 'show', 'partition', partition]
        print(f"Command: {' '.join(scontrol_cmd)}")
        result = subprocess.run(scontrol_cmd, capture_output=True, text=True, timeout=10)
        
        if result.returncode == 0:
            print(result.stdout)
            
            # Look for TRESPerNode or other GPU info
            for line in result.stdout.split('\n'):
                if 'tres' in line.lower() or 'gpu' in line.lower() or 'gres' in line.lower():
                    print(f">>> RELEVANT LINE: {line}")
    except Exception as e:
        print(f"scontrol failed: {e}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python diagnose_sinfo.py <partition_name>")
        print("\nExample:")
        print("  python diagnose_sinfo.py boost_usr_prod")
        sys.exit(1)
    
    partition = sys.argv[1]
    run_sinfo_variants(partition)
    
    print(f"\n{'='*80}")
    print("NEXT STEPS")
    print(f"{'='*80}")
    print("1. Find which command shows GPU information correctly")
    print("2. Note the exact format of the GRES/GPU field")
    print("3. Update the framework's sinfo query to use that format")
    print("4. Update the GRES parser to handle that format")
    print(f"{'='*80}\n")
