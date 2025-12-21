#!/usr/bin/env python3
"""
MPI Runtime Detection and Capability Checking

Detects MPI implementation and version, then determines which flags are supported.
Supports: Intel MPI, OpenMPI, MPICH, Cray MPI, IBM Spectrum MPI

Key issue: Intel MPI does NOT support --report-bindings flag
"""

import subprocess
import re
import logging

logger = logging.getLogger(__name__)


class MPIImplementation(object):
    """MPI implementation information."""
    
    def __init__(self, vendor, version, variant=None):
        self.vendor = vendor  # 'intel', 'openmpi', 'mpich', 'cray', 'spectrum'
        self.version = version  # e.g., '2021.12.1'
        self.variant = variant  # e.g., 'oneapi', 'hpcx'
        
        # Compatibility alias
        self.implementation = vendor  # Alias for vendor
        
        # Capability flags
        self.supports_report_bindings = False
        self.supports_map_by = False
        self.supports_bind_to = False
        self.supports_rank_by = False
        self.launcher_name = 'mpirun'  # or 'mpiexec'
    
    def supports_ppr_mapping(self):
        """Check if this MPI implementation supports ppr (process-per-resource) mapping."""
        return self.supports_map_by
    
    def supports_report_bindings_flag(self):
        """Check if this MPI implementation supports --report-bindings flag."""
        return self.supports_report_bindings
    
    def get_launcher_name(self):
        """Get the launcher command name (mpirun, mpiexec, srun)."""
        return self.launcher_name
    
    def __repr__(self):
        return "MPIImplementation(vendor='{}', version='{}', variant='{}')".format(
            self.vendor, self.version, self.variant or 'default'
        )


class MPIDetector(object):
    """
    Detect MPI implementation and determine supported flags.
    
    Tests actual MPI runtime to determine capabilities.
    """
    
    # Known MPI implementations and their characteristics
    MPI_VENDORS = {
        'intel': {
            'patterns': [r'Intel.*MPI', r'IMPI', r'I_MPI'],
            'supports_report_bindings': False,  # Intel MPI doesn't support this
            'supports_map_by': True,
            'launcher': 'mpirun',
        },
        'openmpi': {
            'patterns': [r'Open\s*MPI', r'OMPI'],
            'supports_report_bindings': True,   # OpenMPI supports this
            'supports_map_by': True,
            'launcher': 'mpirun',
        },
        'mpich': {
            'patterns': [r'MPICH', r'MVAPICH'],
            'supports_report_bindings': False,
            'supports_map_by': False,
            'launcher': 'mpiexec',
        },
        'cray': {
            'patterns': [r'Cray\s*MPI'],
            'supports_report_bindings': False,
            'supports_map_by': False,
            'launcher': 'srun',
        },
        'spectrum': {
            'patterns': [r'IBM\s*Spectrum\s*MPI', r'Platform\s*MPI'],
            'supports_report_bindings': False,
            'supports_map_by': True,
            'launcher': 'mpirun',
        },
    }
    
    def __init__(self):
        self.detected_mpi = None
        self._detection_attempted = False
    
    def detect(self):
        """
        Detect MPI implementation.
        
        Returns:
            MPIImplementation object or None
        """
        self._detection_attempted = True
        
        # Try multiple detection methods
        mpi_info = None
        
        # Method 1: mpirun --version
        mpi_info = self._detect_from_command('mpirun', '--version')
        if mpi_info:
            self.detected_mpi = mpi_info
            logger.info("✓ Detected MPI: {} {} {}".format(
                mpi_info.vendor, mpi_info.version,
                '({})'.format(mpi_info.variant) if mpi_info.variant else ''
            ))
            return mpi_info
        
        # Method 2: mpiexec --version
        mpi_info = self._detect_from_command('mpiexec', '--version')
        if mpi_info:
            self.detected_mpi = mpi_info
            logger.info("✓ Detected MPI: {} {}".format(mpi_info.vendor, mpi_info.version))
            return mpi_info
        
        # Method 3: Check environment variables
        mpi_info = self._detect_from_environment()
        if mpi_info:
            self.detected_mpi = mpi_info
            logger.info("✓ Detected MPI from environment: {} {}".format(
                mpi_info.vendor, mpi_info.version
            ))
            return mpi_info
        
        # Method 4: Module system check
        mpi_info = self._detect_from_modules()
        if mpi_info:
            self.detected_mpi = mpi_info
            logger.info("✓ Detected MPI from modules: {} {}".format(
                mpi_info.vendor, mpi_info.version
            ))
            return mpi_info
        
        logger.warning("⚠ Could not detect MPI implementation, using safe defaults")
        return self._get_safe_default()
    
    def _detect_from_command(self, command, flag):
        """Detect MPI from command output."""
        try:
            result = subprocess.run(
                [command, flag],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode != 0:
                return None
            
            output = result.stdout + result.stderr
            return self._parse_mpi_version(output)
            
        except Exception as e:
            logger.debug("MPI detection via {} {} failed: {}".format(command, flag, e))
            return None
    
    def _detect_from_environment(self):
        """Detect MPI from environment variables."""
        import os
        
        # Intel MPI
        if 'I_MPI_ROOT' in os.environ or 'INTEL_MPI_ROOT' in os.environ:
            version = os.environ.get('I_MPI_VERSION', 'unknown')
            return self._create_mpi_info('intel', version)
        
        # OpenMPI
        if 'OMPI_VERSION' in os.environ:
            version = os.environ.get('OMPI_VERSION', 'unknown')
            return self._create_mpi_info('openmpi', version)
        
        # MPICH
        if 'MPICH_VERSION' in os.environ:
            version = os.environ.get('MPICH_VERSION', 'unknown')
            return self._create_mpi_info('mpich', version)
        
        return None
    
    def _detect_from_modules(self):
        """Detect MPI from loaded modules."""
        try:
            result = subprocess.run(
                ['module', 'list'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            output = result.stdout + result.stderr
            
            # Check for Intel MPI
            if re.search(r'intel.*mpi|impi', output, re.IGNORECASE):
                version_match = re.search(r'(\d+\.\d+\.\d+)', output)
                version = version_match.group(1) if version_match else 'unknown'
                return self._create_mpi_info('intel', version)
            
            # Check for OpenMPI
            if re.search(r'openmpi', output, re.IGNORECASE):
                version_match = re.search(r'(\d+\.\d+\.\d+)', output)
                version = version_match.group(1) if version_match else 'unknown'
                return self._create_mpi_info('openmpi', version)
            
        except Exception as e:
            logger.debug("Module list check failed: {}".format(e))
        
        return None
    
    def _parse_mpi_version(self, output):
        """Parse MPI vendor and version from output."""
        # Check each vendor
        for vendor, info in self.MPI_VENDORS.items():
            for pattern in info['patterns']:
                if re.search(pattern, output, re.IGNORECASE):
                    # Extract version
                    version_patterns = [
                        r'Version\s+(\d+\.\d+\.\d+)',
                        r'(\d+\.\d+\.\d+)',
                    ]
                    
                    version = 'unknown'
                    for vpat in version_patterns:
                        match = re.search(vpat, output)
                        if match:
                            version = match.group(1)
                            break
                    
                    return self._create_mpi_info(vendor, version)
        
        return None
    
    def _create_mpi_info(self, vendor, version, variant=None):
        """Create MPIImplementation object with capabilities."""
        if vendor not in self.MPI_VENDORS:
            return None
        
        info = MPIImplementation(vendor, version, variant)
        vendor_info = self.MPI_VENDORS[vendor]
        
        info.supports_report_bindings = vendor_info['supports_report_bindings']
        info.supports_map_by = vendor_info['supports_map_by']
        info.launcher_name = vendor_info['launcher']
        
        return info
    
    def _get_safe_default(self):
        """Return safe default MPI configuration."""
        # Use conservative settings that work with most MPI implementations
        info = MPIImplementation('generic', 'unknown')
        info.supports_report_bindings = False  # Conservative
        info.supports_map_by = True
        info.launcher_name = 'mpirun'
        return info
    
    def test_flag_support(self, flag):
        """
        Test if a specific flag is supported by running a test command.
        
        Args:
            flag: Flag to test (e.g., '--report-bindings')
        
        Returns:
            True if supported, False otherwise
        """
        try:
            # Try with help or a no-op command
            result = subprocess.run(
                ['mpirun', flag, '--help'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # If it doesn't complain about unrecognized argument, it's likely supported
            output = result.stdout + result.stderr
            if 'unrecognized' in output.lower() or 'invalid' in output.lower():
                return False
            
            return True
            
        except Exception:
            return False
    
    def get_safe_mpirun_flags(self):
        """
        Get safe mpirun flags based on detected MPI.
        
        Returns:
            dict with supported flags
        """
        if not self.detected_mpi:
            self.detect()
        
        mpi = self.detected_mpi or self._get_safe_default()
        
        flags = {
            'report_bindings': '--report-bindings' if mpi.supports_report_bindings else None,
            'map_by_supported': mpi.supports_map_by,
            'launcher': mpi.launcher_name,
        }
        
        logger.debug("Safe MPI flags: {}".format(flags))
        return flags


def test_mpi_detector():
    """Test the MPI detector."""
    import sys
    
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    print("="*70)
    print("MPI Implementation Detector - Test")
    print("="*70)
    print()
    
    detector = MPIDetector()
    mpi = detector.detect()
    
    if mpi:
        print("Detected MPI Implementation:")
        print("  Vendor:  {}".format(mpi.vendor))
        print("  Version: {}".format(mpi.version))
        if mpi.variant:
            print("  Variant: {}".format(mpi.variant))
        print()
        
        print("Supported Features:")
        print("  --report-bindings: {}".format('YES ✓' if mpi.supports_report_bindings else 'NO ✗'))
        print("  --map-by:          {}".format('YES ✓' if mpi.supports_map_by else 'NO ✗'))
        print("  Launcher:          {}".format(mpi.launcher_name))
        print()
        
        print("Safe mpirun command flags:")
        flags = detector.get_safe_mpirun_flags()
        if flags['report_bindings']:
            print("  ✓ Can use: {}".format(flags['report_bindings']))
        else:
            print("  ✗ Cannot use: --report-bindings (will be omitted from commands)")
        
        if flags['map_by_supported']:
            print("  ✓ Can use: --map-by")
        else:
            print("  ✗ Cannot use: --map-by")
        
        print()
        print("="*70)
        print("Recommended Command Format:")
        print("="*70)
        
        if mpi.vendor == 'intel':
            print("For Intel MPI:")
            print("  mpirun -np 4 -ppn 4 $BINARY/app")
            print()
            print("  ⚠ Intel MPI does NOT support --report-bindings")
            print("  ⚠ This flag will cause job failure!")
        elif mpi.vendor == 'openmpi':
            print("For OpenMPI:")
            print("  mpirun -np 4 --map-by ppr:4:node --bind-to core --cpus-per-proc 8 --report-bindings $BINARY/app")
        else:
            print("For {} MPI:".format(mpi.vendor))
            print("  Use conservative command without optional flags")
        
    else:
        print("Could not detect MPI implementation")
        print("Using safe defaults (no --report-bindings)")
    
    print()


if __name__ == '__main__':
    test_mpi_detector()
