"""
Tests for configuration parsing and validation.
"""

import pytest
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.config_parser import load_yaml_config


def test_config_import():
    """Test that config module can be imported."""
    from core import config
    assert config is not None


def test_basic_yaml_structure():
    """Test parsing of basic YAML structure."""
    yaml_content = """
    repository: /path/to/code
    scaling: strong
    nodes: 4
    hardware: cpu
    procs_per_node: 128
    partition: compute
    account: test_account
    time_limit: "01:00:00"
    """
    
    # Write temporary config
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        temp_path = f.name
    
    try:
        config = load_yaml_config(Path(temp_path))
        assert config is not None
        assert config.get('repository') == '/path/to/code'
        assert config.get('scaling') == 'strong'
        assert config.get('nodes') == 4
    finally:
        Path(temp_path).unlink()


def test_scaling_types():
    """Test scaling type validation."""
    valid_types = ['strong', 'weak']
    
    for scaling_type in valid_types:
        yaml_content = f"""
        repository: /path/to/code
        scaling: {scaling_type}
        nodes: 4
        """
        
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            temp_path = f.name
        
        try:
            config = load_yaml_config(Path(temp_path))
            assert config.get('scaling') == scaling_type
        finally:
            Path(temp_path).unlink()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
