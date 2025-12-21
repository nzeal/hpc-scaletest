"""
Test suite for GPU enhancement modules.

Tests source acquisition, build discovery, and binary detection.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import sys
import os

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.source_manager import SourceManager, SourceMetadata
from utils.build_discovery import BuildDiscovery, BuildInfo
from utils.binary_detector import BinaryDetector, BinaryCandidate


class TestSourceManager:
    """Tests for SourceManager."""

    def test_is_git_url_https(self):
        """Test detection of HTTPS Git URLs."""
        manager = SourceManager()
        assert manager._is_git_url("https://github.com/user/repo.git")
        assert manager._is_git_url("https://github.com/user/repo")
        assert manager._is_git_url("https://gitlab.com/user/repo")

    def test_is_git_url_ssh(self):
        """Test detection of SSH Git URLs."""
        manager = SourceManager()
        assert manager._is_git_url("git@github.com:user/repo.git")
        assert manager._is_git_url("git@gitlab.com:user/repo.git")

    def test_is_git_url_local_path(self):
        """Test that local paths are not detected as Git URLs."""
        manager = SourceManager()
        assert not manager._is_git_url("/home/user/repo")
        assert not manager._is_git_url("./repo")
        assert not manager._is_git_url("../repo")

    def test_extract_project_name(self):
        """Test extraction of project name from Git URL."""
        manager = SourceManager()
        assert manager._extract_project_name("https://github.com/user/iPIC3D.git") == "iPIC3D"
        assert manager._extract_project_name("git@github.com:user/gromacs") == "gromacs"
        assert manager._extract_project_name("https://gitlab.com/org/my-app.git") == "my_app"

    def test_acquire_from_local_valid_path(self):
        """Test acquiring source from valid local path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir) / "test_source"
            source_path.mkdir()
            (source_path / "README.md").touch()

            manager = SourceManager()
            result = manager._acquire_from_local(str(source_path))
            assert result == source_path.resolve()

    def test_acquire_from_local_invalid_path(self):
        """Test acquiring source from invalid local path."""
        manager = SourceManager()
        with pytest.raises(ValueError):
            manager._acquire_from_local("/nonexistent/path")

    def test_validate_source_valid(self):
        """Test validation of valid source code."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir)
            (source_path / "CMakeLists.txt").touch()
            (source_path / "README.md").touch()

            manager = SourceManager()
            is_valid, message = manager.validate_source(source_path)
            assert is_valid
            assert "Valid" in message

    def test_validate_source_invalid(self):
        """Test validation of invalid source code."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir)

            manager = SourceManager()
            is_valid, message = manager.validate_source(source_path)
            assert not is_valid
            assert "No build system" in message or "not found" in message

    def test_get_source_metadata(self):
        """Test extraction of source metadata."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir) / "test_project"
            source_path.mkdir()
            (source_path / "README.md").touch()
            (source_path / "CMakeLists.txt").touch()

            manager = SourceManager()
            metadata = manager.get_source_metadata(source_path)

            assert metadata.project_name == "test_project"
            assert metadata.source_path == source_path
            assert metadata.has_readme
            assert metadata.has_cmake
            assert not metadata.is_git_repo


class TestBuildDiscovery:
    """Tests for BuildDiscovery."""

    def test_detect_build_system_cmake(self):
        """Test detection of CMake build system."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir)
            (source_path / "CMakeLists.txt").touch()

            discovery = BuildDiscovery()
            build_system = discovery._detect_build_system(source_path)
            assert build_system == "cmake"

    def test_detect_build_system_make(self):
        """Test detection of Make build system."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir)
            (source_path / "Makefile").touch()

            discovery = BuildDiscovery()
            build_system = discovery._detect_build_system(source_path)
            assert build_system == "make"

    def test_detect_build_system_autotools(self):
        """Test detection of Autotools build system."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir)
            (source_path / "configure").touch()

            discovery = BuildDiscovery()
            build_system = discovery._detect_build_system(source_path)
            assert build_system == "autotools"

    def test_extract_modules_from_readme(self):
        """Test extraction of modules from README."""
        readme_content = """
# Installation

Load required modules:
```
module load cuda/11.8
module load gcc/11.3
module load openmpi/4.1
```
"""
        discovery = BuildDiscovery()
        modules = discovery._extract_modules(readme_content)

        assert "cuda/11.8" in modules
        assert "gcc/11.3" in modules
        assert "openmpi/4.1" in modules

    def test_extract_compiler_requirements(self):
        """Test extraction of compiler requirements."""
        readme_content = """
# Requirements

- CUDA 11.8
- GCC 11.3
- CMake 3.18+
"""
        discovery = BuildDiscovery()
        requirements = discovery._extract_compiler_requirements(readme_content)

        assert requirements.get("cuda") == "11.8"
        assert requirements.get("gcc") == "11.3"
        assert requirements.get("cmake") == "3.18"

    def test_extract_mpi_requirement(self):
        """Test extraction of MPI requirement."""
        discovery = BuildDiscovery()

        # Test OpenMPI
        assert discovery._extract_mpi_requirement("Uses OpenMPI") == "openmpi"

        # Test MPICH
        assert discovery._extract_mpi_requirement("Requires MPICH") == "mpich"

        # Test generic MPI
        assert discovery._extract_mpi_requirement("Needs MPI support") == "mpi"

        # Test no MPI
        assert discovery._extract_mpi_requirement("No MPI required") is None

    def test_analyze_source_with_readme(self):
        """Test complete source analysis with README."""
        with tempfile.TemporaryDirectory() as tmpdir:
            source_path = Path(tmpdir)

            # Create README with build info
            readme_content = """
# Build Instructions

## Requirements
- CUDA 11.8
- GCC 11.3
- OpenMPI 4.1

## Building
```
mkdir build
cd build
cmake ..
make -j
```
"""
            (source_path / "README.md").write_text(readme_content)
            (source_path / "CMakeLists.txt").touch()

            discovery = BuildDiscovery()
            build_info = discovery.analyze_source(source_path)

            assert build_info.build_system == "cmake"
            assert "cuda/11.8" in build_info.required_modules
            assert "gcc/11.3" in build_info.required_modules
            assert build_info.mpi_requirement == "openmpi"
            assert build_info.confidence > 0.5


class TestBinaryDetector:
    """Tests for BinaryDetector."""

    def test_should_skip_file(self):
        """Test file skipping logic."""
        detector = BinaryDetector()

        # Should skip
        assert detector._should_skip_file(Path("test.o"))
        assert detector._should_skip_file(Path("lib.a"))
        assert detector._should_skip_file(Path("lib.so"))
        assert detector._should_skip_file(Path("source.cpp"))
        assert detector._should_skip_file(Path("CMakeFiles/test"))

        # Should not skip
        assert not detector._should_skip_file(Path("iPIC3D"))
        assert not detector._should_skip_file(Path("gromacs"))
        assert not detector._should_skip_file(Path("app"))

    def test_is_elf_executable_valid(self):
        """Test ELF header validation with valid binary."""
        detector = BinaryDetector()

        with tempfile.TemporaryDirectory() as tmpdir:
            binary_path = Path(tmpdir) / "test_binary"

            # Create a minimal valid ELF header
            # ELF magic number + minimal header
            elf_header = b'\x7fELF'  # Magic number
            elf_header += b'\x02'     # 64-bit
            elf_header += b'\x01'     # Little-endian
            elf_header += b'\x01'     # ELF version
            elf_header += b'\x00' * 9  # Padding
            elf_header += b'\x02\x00'  # e_type = executable (2)
            elf_header += b'\x00' * 2  # Padding

            binary_path.write_bytes(elf_header)

            assert detector._is_elf_executable(binary_path)

    def test_is_elf_executable_invalid(self):
        """Test ELF header validation with invalid file."""
        detector = BinaryDetector()

        with tempfile.TemporaryDirectory() as tmpdir:
            binary_path = Path(tmpdir) / "not_elf"
            binary_path.write_text("This is not an ELF file")

            assert not detector._is_elf_executable(binary_path)

    def test_score_location(self):
        """Test location scoring."""
        detector = BinaryDetector()
        source_dir = Path("/home/user/project")

        # Test build/bin location
        score = detector._score_location(Path("/home/user/project/build/bin/app"), source_dir)
        assert score == 25.0

        # Test build location
        score = detector._score_location(Path("/home/user/project/build/app"), source_dir)
        assert score == 20.0

        # Test root location
        score = detector._score_location(Path("/home/user/project/app"), source_dir)
        assert score == 0.0

    def test_score_name(self):
        """Test name scoring."""
        detector = BinaryDetector()
        source_name = "ipic3d"

        # Exact match
        score = detector._score_name(Path("ipic3d"), source_name)
        assert score == 40.0

        # Partial match
        score = detector._score_name(Path("iPIC3D_gpu"), source_name)
        assert score >= 30.0

        # Common pattern match
        score = detector._score_name(Path("gromacs"), "myapp")
        assert score >= 20.0

    def test_score_size(self):
        """Test size scoring."""
        detector = BinaryDetector()

        # Large binary (> 1MB)
        score = detector._score_size(2 * 1024 * 1024)
        assert score == 10.0

        # Medium binary (100KB - 1MB)
        score = detector._score_size(500 * 1024)
        assert score == 5.0

        # Small binary (< 100KB)
        score = detector._score_size(50 * 1024)
        assert score == 0.0

    def test_detect_binary_with_candidates(self):
        """Test binary detection with multiple candidates."""
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "build"
            build_dir.mkdir()
            (build_dir / "bin").mkdir()

            source_dir = Path(tmpdir) / "source"
            source_dir.mkdir()

            # Create test binaries
            best_binary = build_dir / "bin" / "iPIC3D"
            best_binary.write_bytes(b'\x7fELF' + b'\x00' * 16)
            best_binary.chmod(0o755)

            other_binary = build_dir / "test"
            other_binary.write_bytes(b'\x7fELF' + b'\x00' * 16)
            other_binary.chmod(0o755)

            detector = BinaryDetector()
            detected = detector.detect_binary(build_dir, source_dir)

            # Should prefer binary in bin directory
            assert detected is not None
            assert "bin" in str(detected)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
