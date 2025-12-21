"""
GPU Allocation Directives

Provides scheduler-agnostic GPU allocation directives.
Automatically translates to scheduler-specific syntax.
"""

import logging
from core.decorators import register_directive

logger = logging.getLogger(__name__)


@register_directive("gpu")
def gpu_directive(num_gpus, gpu_type=None, gpu_bind=None, **kwargs):
    """
    Generate GPU allocation directive for scheduler.
    
    Parameters:
        num_gpus (int): Number of GPUs to allocate per node
        gpu_type (str, optional): Specific GPU type (e.g., 'a100', 'v100', 'mi250')
        gpu_bind (str, optional): GPU binding strategy ('closest', 'single', 'none')
        **kwargs: Additional scheduler-specific options
    
    Returns:
        dict: Scheduler-specific directives
            {
                'slurm': '--gpus-per-node=4 --gpu-bind=closest',
                'pbs': '-l select=1:ncpus=32:ngpus=4',
                'lsf': 'bsub -gpu "num=4:mode=shared"'
            }
    
    Examples:
        >>> gpu_directive(4)
        {'slurm': '--gpus-per-node=4', ...}
        
        >>> gpu_directive(4, gpu_type='a100', gpu_bind='closest')
        {'slurm': '--gpus-per-node=a100:4 --gpu-bind=closest', ...}
    """
    directives = {
        'slurm': _slurm_gpu(num_gpus, gpu_type, gpu_bind, **kwargs),
        'pbs': _pbs_gpu(num_gpus, gpu_type, gpu_bind, **kwargs),
        'lsf': _lsf_gpu(num_gpus, gpu_type, gpu_bind, **kwargs),
        'sge': _sge_gpu(num_gpus, gpu_type, gpu_bind, **kwargs),
    }
    
    logger.debug(f"Generated GPU directives for {num_gpus} GPUs "
                f"(type={gpu_type}, bind={gpu_bind})")
    
    return directives


def _slurm_gpu(num_gpus, gpu_type, gpu_bind, **kwargs):
    """
    SLURM-specific GPU directive.
    
    SLURM GPU options:
        --gpus=N                  # Total GPUs across all nodes
        --gpus-per-node=N         # GPUs per node
        --gpus-per-task=N         # GPUs per MPI task
        --gres=gpu:type:N         # Legacy format
        --gpu-bind=closest        # Bind to closest GPU
    """
    parts = []
    
    # GPU allocation
    if gpu_type:
        # Specific GPU type: --gpus-per-node=a100:4
        parts.append(f"--gpus-per-node={gpu_type}:{num_gpus}")
    else:
        # Any GPU: --gpus-per-node=4
        parts.append(f"--gpus-per-node={num_gpus}")
    
    # GPU binding
    if gpu_bind:
        parts.append(f"--gpu-bind={gpu_bind}")
    
    # Additional SLURM-specific options
    if kwargs.get('gpu_freq'):
        parts.append(f"--gpu-freq={kwargs['gpu_freq']}")
    
    if kwargs.get('gpus_per_task'):
        parts.append(f"--gpus-per-task={kwargs['gpus_per_task']}")
    
    return ' '.join(parts)


def _pbs_gpu(num_gpus, gpu_type, gpu_bind, **kwargs):
    """
    PBS/Torque-specific GPU directive.
    
    PBS GPU options:
        -l select=1:ncpus=32:ngpus=4
        -l select=1:ncpus=32:ngpus=4:gpu_model=a100
    """
    # PBS uses resource lists
    resource_parts = [f"ngpus={num_gpus}"]
    
    if gpu_type:
        resource_parts.append(f"gpu_model={gpu_type}")
    
    # Note: PBS doesn't have direct gpu-bind, relies on system placement
    
    return f"-l select=1:{':'.join(resource_parts)}"


def _lsf_gpu(num_gpus, gpu_type, gpu_bind, **kwargs):
    """
    LSF-specific GPU directive.
    
    LSF GPU options:
        -gpu "num=4"
        -gpu "num=4:mode=shared:mps=yes"
        -R "select[gpu_model==a100]"
    """
    gpu_parts = [f"num={num_gpus}"]
    
    # GPU sharing mode
    mode = kwargs.get('gpu_mode', 'shared')  # 'shared' or 'exclusive_process'
    gpu_parts.append(f"mode={mode}")
    
    # MPS support
    if kwargs.get('mps', False):
        gpu_parts.append("mps=yes")
    
    gpu_directive = f'-gpu "{":".join(gpu_parts)}"'
    
    # GPU type selection
    if gpu_type:
        gpu_directive += f' -R "select[gpu_model=={gpu_type}]"'
    
    return gpu_directive


def _sge_gpu(num_gpus, gpu_type, gpu_bind, **kwargs):
    """
    SGE (Sun Grid Engine) GPU directive.
    
    SGE GPU options:
        -l gpu=4
        -l gpu_type=a100
    """
    parts = [f"-l gpu={num_gpus}"]
    
    if gpu_type:
        parts.append(f"-l gpu_type={gpu_type}")
    
    return ' '.join(parts)


# Additional helper functions

def validate_gpu_allocation(num_gpus, available_gpus):
    """
    Validate GPU allocation request.
    
    Args:
        num_gpus: Requested GPUs
        available_gpus: Available GPUs per node
    
    Returns:
        tuple: (is_valid, error_message)
    """
    if num_gpus <= 0:
        return False, "Number of GPUs must be positive"
    
    if num_gpus > available_gpus:
        return False, f"Requested {num_gpus} GPUs but only {available_gpus} available"
    
    return True, ""


def get_gpu_binding_strategies():
    """
    Get available GPU binding strategies.
    
    Returns:
        dict: Mapping of strategy name to description
    """
    return {
        'closest': 'Bind to closest GPU based on CPU affinity',
        'single': 'Bind to single GPU (exclusive)',
        'multiple': 'Bind to multiple GPUs (shared)',
        'none': 'No explicit binding',
        'map_gpu': 'Custom GPU mapping',
    }


# Register for auto-discovery
__all__ = ['gpu_directive', 'validate_gpu_allocation', 'get_gpu_binding_strategies']
