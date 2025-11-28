"""
Common parsing utilities.

Consolidates duplicate parsing logic from multiple modules.
"""

from typing import Any, Tuple, List, Union


def parse_tuple(value: Any) -> Tuple:
    """
    Parse a tuple value from various formats.
    
    Supports:
    - Lists: [1, 2, 3]
    - Tuples: (1, 2, 3)
    - Comma-separated strings: "1,2,3"
    - Space-separated strings: "1 2 3"
    
    Args:
        value: Input value in various formats
        
    Returns:
        Tuple of parsed values
        
    Examples:
        >>> parse_tuple([1, 2, 3])
        (1, 2, 3)
        >>> parse_tuple("1.5, 2.0, 3.5")
        (1.5, 2.0, 3.5)
        >>> parse_tuple("14 8 1")
        (14, 8, 1)
    """
    if isinstance(value, (list, tuple)):
        return tuple(value)
    
    elif isinstance(value, str):
        # Try comma-separated first
        if ',' in value:
            parts = [x.strip() for x in value.split(',')]
        else:
            # Space-separated
            parts = value.split()
        
        # Convert to appropriate types
        result = []
        for part in parts:
            try:
                # Try int first
                if '.' not in part and 'e' not in part.lower():
                    result.append(int(part))
                else:
                    # Float
                    result.append(float(part))
            except ValueError:
                # Keep as string
                result.append(part)
        
        return tuple(result)
    
    elif isinstance(value, dict):
        # Dictionary format: {x: 1, y: 2, z: 3}
        return tuple(value.get(k) for k in ['x', 'y', 'z'] if k in value)
    
    else:
        raise ValueError(f"Cannot parse tuple from: {value} (type: {type(value)})")


def parse_list(value: Any) -> List:
    """
    Parse a list value from various formats.
    
    Args:
        value: Input value
        
    Returns:
        List of parsed values
    """
    if isinstance(value, list):
        return value
    elif isinstance(value, tuple):
        return list(value)
    elif isinstance(value, str):
        return list(parse_tuple(value))
    else:
        raise ValueError(f"Cannot parse list from: {value}")


def parse_numeric(value: Any) -> Union[int, float]:
    """
    Parse a numeric value.
    
    Args:
        value: Input value
        
    Returns:
        Parsed number (int or float)
    """
    if isinstance(value, (int, float)):
        return value
    
    elif isinstance(value, str):
        value = value.strip()
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            else:
                return int(value)
        except ValueError:
            raise ValueError(f"Cannot parse numeric from: {value}")
    
    else:
        raise ValueError(f"Cannot parse numeric from: {value} (type: {type(value)})")


def parse_bool(value: Any) -> bool:
    """
    Parse a boolean value.
    
    Handles: true/false, yes/no, 1/0, True/False
    
    Args:
        value: Input value
        
    Returns:
        Boolean value
    """
    if isinstance(value, bool):
        return value
    
    elif isinstance(value, str):
        value = value.lower().strip()
        if value in ['true', 'yes', '1', 'on', 'enabled']:
            return True
        elif value in ['false', 'no', '0', 'off', 'disabled']:
            return False
        else:
            raise ValueError(f"Cannot parse bool from: {value}")
    
    elif isinstance(value, (int, float)):
        return bool(value)
    
    else:
        raise ValueError(f"Cannot parse bool from: {value}")
