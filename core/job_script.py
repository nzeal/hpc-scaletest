#!/usr/bin/env python3
"""
Job Script Generator Module

Generates SLURM job scripts with correct resource allocation and MPI launch
commands based on automatically detected hardware topology.

Design Principles:
1. NO HARDCODED VALUES - All configuration derived from topology
2. CORRECT RESOURCE MAPPING - SLURM allocates resources, MPI uses them
3. AUTOMATIC GPU BINDING - Via CUDA_VISIBLE_DEVICES and local rank

This module integrates:
- core/topology.py for hardware detection
- core/mpi_command.py for MPI command generation

Author: HPC-ScaleTest Contributors
"""

import os
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field

from core.topology import (
    NodeTopology, MPIMapping, GPUVendor,
    TopologyDetector, get_topology_detector
)
from core.mpi_command import (
    MPICommandGenerator, MPIDetector, MPIInfo,
    generate_gpu_binding_script
)

logger = logging.getLogger(__name__)


@dataclass
class JobConfiguration:
    """
    Job configuration parameters.
    
    This class holds user-specified job parameters that cannot be
    auto-detected from topology.
    """
    job_name: str = "hpc_scaletest"
    num_nodes: int = 1
    time_limit: str = "01:00:00"
    partition: Optional[str] = None
    account: Optional[str] = None
    qos: Optional[str] = None
    
    # Executable configuration
    executable: str = ""
    executable_args: List[str] = field(default_factory=list)
    
    # Optional overrides (if None, auto-detected)
    ranks_per_node: Optional[int] = None
    cores_per_rank: Optional[int] = None
    gpus_per_node: Optional[int] = None
    
    # Module system
    modules: List[str] = field(default_factory=list)
    
    # Additional SLURM options
    extra_slurm_options: Dict[str, str] = field(default_factory=dict)
    
    # Output configuration
    output_dir: str = "."
    stdout_file: str = "job_%j.out"
    stderr_file: str = "job_%j.err"
    
    # Feature flags
    exclusive: bool = True
    verbose_mpi: bool = False


class JobScriptGenerator:
    """
    Generate SLURM job scripts with automatic topology detection.
    
    This class produces job scripts that correctly allocate SLURM resources
    and launch MPI jobs based on the detected hardware topology.
    
    Key Responsibilities:
    1. Detect hardware topology for the target partition
    2. Compute optimal MPI mapping
    3. Generate SLURM directives
    4. Generate MPI launch command
    5. Generate GPU binding script (if needed)
    
    Usage:
        generator = JobScriptGenerator()
        script = generator.generate(config)
        generator.write(script, "job.sh")
    """
    
    def __init__(self):
        """Initialize the job script generator."""
        self.topology_detector = get_topology_detector()
        self.mpi_detector = MPIDetector()
    
    def generate(self, config: JobConfiguration) -> str:
        """
        Generate a complete SLURM job script.
        
        Args:
            config: Job configuration parameters
        
        Returns:
            Complete job script as a string
        """
        # Step 1: Detect topology
        topology = self.topology_detector.detect(config.partition)
        
        # Apply user overrides
        if config.gpus_per_node is not None:
            # User override - adjust topology
            topology = NodeTopology(
                cpu_cores=topology.cpu_cores,
                cpu_sockets=topology.cpu_sockets,
                cores_per_socket=topology.cores_per_socket,
                threads_per_core=topology.threads_per_core,
                physical_cores=topology.physical_cores,
                gpus=config.gpus_per_node,
                gpu_vendor=topology.gpu_vendor,
                gpu_model=topology.gpu_model,
                memory_gb=topology.memory_gb,
                partition=topology.partition,
                detection_method=topology.detection_method + " + user override",
            )
        
        # Step 2: Compute MPI mapping
        mapping = self.topology_detector.compute_mpi_mapping(
            topology=topology,
            num_nodes=config.num_nodes,
            user_ranks_per_node=config.ranks_per_node,
            user_cores_per_rank=config.cores_per_rank,
        )
        
        # Step 3: Detect MPI implementation
        mpi_info = self.mpi_detector.detect()
        
        # Step 4: Generate script sections
        script_parts = [
            self._generate_header(),
            self._generate_slurm_directives(config, topology, mapping),
            self._generate_environment_section(config, topology),
            self._generate_info_section(topology, mapping),
            self._generate_gpu_binding_section(topology, config.output_dir),
            self._generate_execution_section(config, topology, mapping, mpi_info),
        ]
        
        return '\n'.join(script_parts)
    
    def _generate_header(self) -> str:
        """Generate script header."""
        return '''#!/bin/bash
# =============================================================================
# SLURM Job Script
# Generated by HPC-ScaleTest
#
# Hardware topology and MPI mapping are automatically detected.
# No site-specific constants are hardcoded in this script.
# =============================================================================
'''
    
    def _generate_slurm_directives(
        self,
        config: JobConfiguration,
        topology: NodeTopology,
        mapping: MPIMapping,
    ) -> str:
        """
        Generate SLURM directives based on topology and mapping.
        
        Key insight: SLURM allocates resources, MPI uses them differently.
        
        For GPU jobs:
        - --ntasks-per-node should equal the number of MPI ranks
        - --cpus-per-task should equal cores per rank
        - --gres=gpu:N requests N GPUs
        
        The total CPUs allocated = ntasks-per-node × cpus-per-task
        This should equal the cores we want to use (not necessarily all cores).
        """
        lines = ['# SLURM Directives']
        
        # Job identification
        lines.append(f'#SBATCH --job-name={config.job_name}')
        
        # Resource allocation
        lines.append(f'#SBATCH --nodes={config.num_nodes}')
        lines.append(f'#SBATCH --ntasks-per-node={mapping.ranks_per_node}')
        lines.append(f'#SBATCH --cpus-per-task={mapping.cores_per_rank}')
        
        # GPU allocation (if GPUs present)
        if topology.gpus > 0:
            lines.append(f'#SBATCH --gres=gpu:{topology.gpus}')
        
        # Time limit
        lines.append(f'#SBATCH --time={config.time_limit}')
        
        # Partition (required if detected or specified)
        if config.partition:
            lines.append(f'#SBATCH --partition={config.partition}')
        elif topology.partition:
            lines.append(f'#SBATCH --partition={topology.partition}')
        
        # Account (optional)
        if config.account:
            lines.append(f'#SBATCH --account={config.account}')
        
        # QOS (optional)
        if config.qos:
            lines.append(f'#SBATCH --qos={config.qos}')
        
        # Exclusive mode
        if config.exclusive:
            lines.append('#SBATCH --exclusive')
        
        # Output files
        lines.append(f'#SBATCH --output={config.stdout_file}')
        lines.append(f'#SBATCH --error={config.stderr_file}')
        
        # Extra options
        for key, value in config.extra_slurm_options.items():
            lines.append(f'#SBATCH --{key}={value}')
        
        lines.append('')
        return '\n'.join(lines)
    
    def _generate_environment_section(
        self,
        config: JobConfiguration,
        topology: NodeTopology,
    ) -> str:
        """Generate environment setup section."""
        lines = ['# Environment Setup']
        
        # Set common environment variables
        lines.append('set -e  # Exit on error')
        lines.append('')
        
        # Module loading
        if config.modules:
            lines.append('# Load required modules')
            for module in config.modules:
                lines.append(f'module load {module}')
            lines.append('')
        
        # GPU-specific environment
        if topology.gpus > 0:
            lines.append('# GPU Environment')
            if topology.gpu_vendor == GPUVendor.NVIDIA:
                lines.append('# CUDA_VISIBLE_DEVICES will be set by bind_gpu.sh')
            elif topology.gpu_vendor == GPUVendor.AMD:
                lines.append('# ROCR_VISIBLE_DEVICES will be set by bind_gpu.sh')
            elif topology.gpu_vendor == GPUVendor.INTEL:
                lines.append('# ZE_AFFINITY_MASK will be set by bind_gpu.sh')
            lines.append('')
        
        return '\n'.join(lines)
    
    def _generate_info_section(
        self,
        topology: NodeTopology,
        mapping: MPIMapping,
    ) -> str:
        """Generate runtime information section."""
        lines = ['# Job Information']
        lines.append('echo "========================================"')
        lines.append('echo "HPC-ScaleTest Job Starting"')
        lines.append('echo "========================================"')
        lines.append('echo "Date: $(date)"')
        lines.append('echo "Hostname: $(hostname)"')
        lines.append('echo "Job ID: $SLURM_JOB_ID"')
        lines.append(f'echo "Nodes: $SLURM_JOB_NUM_NODES"')
        lines.append('')
        
        # Hardware info
        lines.append('echo ""')
        lines.append('echo "Hardware Configuration (auto-detected):"')
        lines.append(f'echo "  CPU cores per node: {topology.cpu_cores}"')
        lines.append(f'echo "  GPUs per node: {topology.gpus}"')
        if topology.gpu_model:
            lines.append(f'echo "  GPU model: {topology.gpu_model}"')
        lines.append(f'echo "  Detection method: {topology.detection_method}"')
        lines.append('')
        
        # MPI mapping info
        lines.append('echo ""')
        lines.append('echo "MPI Configuration:"')
        lines.append(f'echo "  Ranks per node: {mapping.ranks_per_node}"')
        lines.append(f'echo "  Cores per rank: {mapping.cores_per_rank}"')
        lines.append(f'echo "  GPUs per rank: {mapping.gpus_per_rank}"')
        lines.append(f'echo "  Total ranks: {mapping.total_ranks}"')
        lines.append('')
        
        # SLURM environment
        lines.append('echo ""')
        lines.append('echo "SLURM Environment:"')
        lines.append('echo "  SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"')
        lines.append('echo "  SLURM_TASKS_PER_NODE: $SLURM_TASKS_PER_NODE"')
        if topology.gpus > 0:
            lines.append('echo "  SLURM_GPUS_ON_NODE: $SLURM_GPUS_ON_NODE"')
            lines.append('echo "  CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"')
        lines.append('echo "========================================"')
        lines.append('echo ""')
        lines.append('')
        
        return '\n'.join(lines)
    
    def _generate_gpu_binding_section(
        self,
        topology: NodeTopology,
        output_dir: str,
    ) -> str:
        """Generate GPU binding script section."""
        if topology.gpus == 0:
            return '# No GPU binding needed (CPU-only job)\n'
        
        # Determine GPU environment variable
        if topology.gpu_vendor == GPUVendor.NVIDIA:
            gpu_env = "CUDA_VISIBLE_DEVICES"
        elif topology.gpu_vendor == GPUVendor.AMD:
            gpu_env = "ROCR_VISIBLE_DEVICES"
        elif topology.gpu_vendor == GPUVendor.INTEL:
            gpu_env = "ZE_AFFINITY_MASK"
        else:
            gpu_env = "CUDA_VISIBLE_DEVICES"
        
        return f'''# GPU Binding Script
# This script ensures each MPI rank gets a unique GPU via OMPI_COMM_WORLD_LOCAL_RANK
cat > bind_gpu.sh << 'EOF'
#!/bin/bash
# Determine local rank
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
else
    LOCAL_RANK=0
fi
export {gpu_env}=$LOCAL_RANK
exec "$@"
EOF
chmod +x bind_gpu.sh

'''
    
    def _generate_execution_section(
        self,
        config: JobConfiguration,
        topology: NodeTopology,
        mapping: MPIMapping,
        mpi_info: MPIInfo,
    ) -> str:
        """Generate the execution section with MPI launch command."""
        lines = ['# Execution']
        lines.append('cd $SLURM_SUBMIT_DIR')
        lines.append('')
        
        # Generate MPI command
        generator = MPICommandGenerator(mpi_info)
        
        # Include GPU binding if GPUs present
        gpu_binding = './bind_gpu.sh' if topology.gpus > 0 else None
        
        cmd = generator.generate(
            topology=topology,
            mapping=mapping,
            executable=config.executable,
            args=config.executable_args,
            num_nodes=config.num_nodes,
            verbose=config.verbose_mpi,
            gpu_binding_script=gpu_binding,
        )
        
        # Format command for script
        cmd_str = ' \\\n    '.join(cmd)
        
        lines.append('echo "Launching MPI job..."')
        lines.append(f'{cmd_str}')
        lines.append('')
        
        # Capture exit code
        lines.append('EXIT_CODE=$?')
        lines.append('')
        
        # Job completion
        lines.append('echo ""')
        lines.append('echo "========================================"')
        lines.append('echo "Job completed at: $(date)"')
        lines.append('echo "Exit code: $EXIT_CODE"')
        lines.append('echo "========================================"')
        lines.append('')
        lines.append('exit $EXIT_CODE')
        
        return '\n'.join(lines)
    
    def write(
        self,
        script: str,
        filepath: str,
        make_executable: bool = True,
    ) -> str:
        """
        Write script to file.
        
        Args:
            script: Script content
            filepath: Output file path
            make_executable: Make the script executable
        
        Returns:
            Absolute path to the written file
        """
        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(path, 'w') as f:
            f.write(script)
        
        if make_executable:
            os.chmod(path, 0o755)
        
        logger.info(f"Generated job script: {path.absolute()}")
        return str(path.absolute())
    
    def generate_and_write(
        self,
        config: JobConfiguration,
        filepath: str,
    ) -> str:
        """
        Generate and write job script in one step.
        
        Args:
            config: Job configuration
            filepath: Output file path
        
        Returns:
            Absolute path to the written file
        """
        script = self.generate(config)
        return self.write(script, filepath)


def generate_job_script(
    executable: str,
    num_nodes: int = 1,
    partition: Optional[str] = None,
    time_limit: str = "01:00:00",
    job_name: str = "hpc_scaletest",
    modules: List[str] = None,
    args: List[str] = None,
    account: Optional[str] = None,
    qos: Optional[str] = None,
    output_file: Optional[str] = None,
    verbose: bool = False,
) -> str:
    """
    Convenience function to generate a job script.
    
    This is the simplest interface for generating a job script.
    Hardware topology is auto-detected based on the partition.
    
    Args:
        executable: Path to the executable
        num_nodes: Number of nodes
        partition: SLURM partition name
        time_limit: Time limit (HH:MM:SS)
        job_name: Job name
        modules: List of modules to load
        args: Executable arguments
        account: SLURM account
        qos: Quality of Service
        output_file: Write script to this file (optional)
        verbose: Enable verbose MPI output
    
    Returns:
        Job script as a string (and writes to file if output_file specified)
    
    Example:
        script = generate_job_script(
            executable="./my_app",
            num_nodes=4,
            partition="gpu",
            modules=["cuda", "openmpi"],
        )
    """
    config = JobConfiguration(
        job_name=job_name,
        num_nodes=num_nodes,
        time_limit=time_limit,
        partition=partition,
        account=account,
        qos=qos,
        executable=executable,
        executable_args=args or [],
        modules=modules or [],
        verbose_mpi=verbose,
    )
    
    generator = JobScriptGenerator()
    script = generator.generate(config)
    
    if output_file:
        generator.write(script, output_file)
    
    return script
