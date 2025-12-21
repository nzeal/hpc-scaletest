"""
Setup script for HPC-ScaleTest package.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

setup(
    name="hpc-scaletest",
    version="2.0.0",
    author="HPC ScaleTest Contributors",
    description="A modular framework for running benchmark scaling tests on HPC systems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/hpc-scaletest",
    packages=find_packages(exclude=["tests", "examples"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: System :: Clustering",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=[
        # Minimal dependencies - most backends use standard library
    ],
    extras_require={
        "yaml": ["pyyaml>=5.0"],
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=22.0",
            "flake8>=5.0",
            "mypy>=0.990",
        ],
    },
    entry_points={
        "console_scripts": [
            "hpc-scaletest=hpc_auto:main",
            "scaletest=hpc_auto:main",  # Backward compatibility
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
