================================================================
Strong Scaling Algorithm: Complete Technical Reference
================================================================

.. contents:: Table of Contents
   :local:
   :depth: 3

Introduction for Everyone
===========================

What is Strong Scaling?
-------------------------

**Strong scaling** is about **making your program run faster** by using more computers.

Think of it like this:

Imagine you need to paint a large house:
- **1 painter**: Takes 10 hours
- **2 painters**: Takes 5 hours (twice as fast!)
- **4 painters**: Takes 2.5 hours (four times as fast!)

The **size of the house doesn't change** - you just use more painters to finish faster.

That's strong scaling: **Same problem size, more workers = faster completion!**

Real-World Example
~~~~~~~~~~~~~~~~~~~

**Scenario**: You have a climate simulation with 1 million grid points.

- On **1 computer**: Takes 24 hours
- On **2 computers**: Takes 12 hours (ideally)
- On **4 computers**: Takes 6 hours (ideally)
- On **8 computers**: Takes 3 hours (ideally)

The simulation size (1 million points) **stays the same**. You're just using more 
computers to finish faster.

HPC-ScaleTest **automatically calculates** how to split this work across any number 
of computers!

Why is Strong Scaling Important?
----------------------------------

**Time is Money**

If your simulation takes 24 hours on 1 computer:
- Using 8 computers → Finish in 3 hours
- Get your results **8 times faster**!
- Run more experiments in less time

**Real Research Examples**

1. **Weather Forecasting**: Must finish before the weather arrives!
2. **Drug Discovery**: Test thousands of molecules quickly
3. **Crash Simulations**: Fast results for design iterations
4. **Genomics**: Process huge datasets in reasonable time

Strong Scaling vs Weak Scaling
--------------------------------

These are the two main ways to use more computers:

.. code-block:: text

   STRONG SCALING (this document):
   ================================
   Problem size: FIXED
   Goal: Make it FASTER
   
   Example:
   - 1 computer: 1M grid points, 10 hours
   - 2 computers: 1M grid points, 5 hours  ← Same problem, half the time
   - 4 computers: 1M grid points, 2.5 hours
   
   
   WEAK SCALING (other document):
   ==============================
   Problem size: GROWS with computers
   Goal: Keep time CONSTANT
   
   Example:
   - 1 computer: 1M grid points, 10 hours
   - 2 computers: 2M grid points, 10 hours  ← Bigger problem, same time
   - 4 computers: 4M grid points, 10 hours

**When to use strong scaling**: You have a fixed problem and want faster results.

**When to use weak scaling**: You want to solve bigger problems without waiting longer.

Conceptual Foundation
======================

How Computers Share Work
--------------------------

When you run a program on multiple computers (called "parallel computing"), 
the computers need to:

1. **Divide the work** - Split the problem into pieces
2. **Work independently** - Each computer handles its piece
3. **Communicate** - Share results with each other
4. **Combine** - Put the final answer together

The Big Challenge: Amdahl's Law
---------------------------------

**Not everything can be parallelized!**

Every program has:
- **Parallel parts**: Can be split across computers (e.g., calculating grid points)
- **Serial parts**: Must be done by one computer (e.g., reading the input file)

**Amdahl's Law** (explained simply):

.. code-block:: text

   If 90% of your program can be parallelized:
   
   - 1 computer:   100 seconds total
   - 2 computers:  55 seconds  (10 serial + 45 parallel)
   - 4 computers:  32.5 seconds (10 serial + 22.5 parallel)
   - 8 computers:  21.25 seconds (10 serial + 11.25 parallel)
   - ∞ computers:  10 seconds (serial part is the limit!)

**Key insight**: The serial part (10%) limits how much speedup you can get!

Efficiency Metrics Explained
------------------------------

When we use more computers, we measure how well it's working:

**Speedup**: How much faster did we get?

.. math::

   \text{Speedup} = \frac{\text{Time on 1 computer}}{\text{Time on N computers}}

**Example**:
- 1 computer: 100 seconds
- 4 computers: 30 seconds
- Speedup = 100 / 30 = 3.33×

**Ideal** would be 4× (4 computers = 4× faster)
**Actual** is 3.33× (pretty good!)

**Efficiency**: How effectively are we using the computers?

.. math::

   \text{Efficiency} = \frac{\text{Speedup}}{\text{Number of computers}} \times 100\%

**Example**:
- Speedup: 3.33×
- Computers: 4
- Efficiency = 3.33 / 4 = 83.25%

**Interpretation**:
- 100% = Perfect! (Each computer contributes fully)
- 80%+ = Good (minor overhead)
- 50-80% = Acceptable (some inefficiency)
- <50% = Poor (too much overhead)

Mathematical Foundations
=========================

Problem Definition
-------------------

Given:
- A computational problem with **fixed size**
- A baseline configuration with :math:`P_0` processes
- A target configuration with :math:`P` processes where :math:`P > P_0`

Goal:
Calculate MPI topology (process grid) that maintains:

1. **Fixed problem size** (grid dimensions don't change)
2. **Proper load distribution** (work evenly divided)
3. **Optimal communication** (minimize processor idle time)

Mathematical Formulation
-------------------------

For a 3D problem with baseline grid :math:`(n_x, n_y, n_z)` on :math:`P_0` processes:

**Baseline Configuration**:

.. math::

   P_0 = p_x^0 \times p_y^0 \times p_z^0

where :math:`(p_x^0, p_y^0, p_z^0)` is the baseline MPI topology.

**Target Configuration**:

.. math::

   P = p_x \times p_y \times p_z

**Constraint**: Grid dimensions remain constant:

.. math::

   n_x, n_y, n_z \text{ are fixed}

**Work per Process**:

.. math::

   W_{\text{baseline}} &= \frac{n_x \times n_y \times n_z}{P_0} \\
   W_{\text{target}} &= \frac{n_x \times n_y \times n_z}{P}

**Key Relationship**:

.. math::

   W_{\text{target}} = W_{\text{baseline}} \times \frac{P_0}{P}

As :math:`P` increases, work per process decreases proportionally (ideal strong scaling).

MPI Topology Calculation
--------------------------

**Goal**: Find :math:`(p_x, p_y, p_z)` such that:

1. :math:`p_x \times p_y \times p_z = P` (correct total)
2. Approximately cubic ratio (optimal communication)
3. Divisibility: :math:`n_x \% p_x = 0`, :math:`n_y \% p_y = 0`, :math:`n_z \% p_z = 0`

**Algorithm**: Prime factorization with balanced distribution

.. math::

   P = 2^a \times 3^b \times 5^c \times \ldots

Distribute prime factors to :math:`(p_x, p_y, p_z)` maintaining:

.. math::

   \frac{p_x}{p_y} \approx \frac{p_y}{p_z} \approx 1

Algorithm Design
=================

Overview
---------

HPC-ScaleTest's strong scaling algorithm has three main steps:

.. code-block:: text

   INPUT: 
   - Baseline configuration (e.g., 8 processes)
   - Target process counts (e.g., 16, 32, 64, 128 processes)
   - Grid dimensions (e.g., 128 × 128 × 128)
   
   ALGORITHM:
   
   Step 1: Calculate MPI Topology
           For each target process count, find best (px, py, pz)
           that multiplies to give the target
   
   Step 2: Verify Divisibility
           Ensure grid can be evenly divided
           (e.g., 128 cells ÷ 4 processes = 32 cells each)
   
   Step 3: Generate Configuration
           Create input files with correct settings
   
   OUTPUT:
   - Configuration files for each process count
   - Job scripts to run experiments
   - Verification that work is evenly distributed

MPI Topology Calculation Algorithm
------------------------------------

**Problem**: Given :math:`P` processes, find :math:`(p_x, p_y, p_z)` where :math:`p_x \times p_y \times p_z = P`

**Method**: Prime Factorization with Balanced Distribution

**Step 1: Factor the number**

.. code-block:: text

   def factorize(n):
       """Find all prime factors of n.
       
       Example: 24 = 2 × 2 × 2 × 3
       Returns: [2, 2, 2, 3]
       """
       factors = []
       d = 2
       while d * d <= n:
           while n % d == 0:
               factors.append(d)
               n //= d
           d += 1
       if n > 1:
           factors.append(n)
       return factors

**Example**: :math:`P = 24`

.. code-block:: text

   24 = 2 × 2 × 2 × 3
   factors = [2, 2, 2, 3]

**Step 2: Distribute factors to create balanced topology**

.. code-block:: text

   def calculate_mpi_topology(num_procs):
       """Calculate balanced 3D MPI topology.
       
       Strategy: Distribute prime factors to maintain cubic ratio.
       Goal: px ≈ py ≈ pz for optimal communication
       """
       if num_procs == 1:
           return (1, 1, 1)
       
       # Get prime factors
       factors = factorize(num_procs)
       
       # Start with (1, 1, 1)
       px, py, pz = 1, 1, 1
       
       # Distribute factors to maintain balance
       for factor in sorted(factors, reverse=True):
           # Find smallest dimension
           if px <= py and px <= pz:
               px *= factor
           elif py <= pz:
               py *= factor
           else:
               pz *= factor
       
       return (px, py, pz)

**Example**: :math:`P = 24`, factors = [3, 2, 2, 2]

.. code-block:: text

   Start: (1, 1, 1)
   
   Factor 3 (largest): 
   - Smallest is px=1, so px = 1 × 3 = 3
   - Result: (3, 1, 1)
   
   Factor 2:
   - Smallest is py=1, so py = 1 × 2 = 2
   - Result: (3, 2, 1)
   
   Factor 2:
   - Smallest is pz=1, so pz = 1 × 2 = 2
   - Result: (3, 2, 2)
   
   Factor 2:
   - Smallest is py=2, so py = 2 × 2 = 4
   - Result: (3, 4, 2)
   
   Final: (3, 4, 2)
   Verify: 3 × 4 × 2 = 24 ✓

Divisibility Verification
---------------------------

**Critical Check**: Can the grid be evenly divided?

.. code-block:: python

   def verify_divisibility(grid_dims, mpi_topo):
       """Verify grid can be evenly divided among processes.
       
       Args:
           grid_dims: (nx, ny, nz) - Total grid size
           mpi_topo: (px, py, pz) - Process grid
       
       Returns:
           True if evenly divisible, False otherwise
       """
       nx, ny, nz = grid_dims
       px, py, pz = mpi_topo
       
       # Check each dimension
       if nx % px != 0:
           print(f"Warning: {nx} cells not divisible by {px} processes in X")
           return False
       if ny % py != 0:
           print(f"Warning: {ny} cells not divisible by {py} processes in Y")
           return False
       if nz % pz != 0:
           print(f"Warning: {nz} cells not divisible by {pz} processes in Z")
           return False
       
       return True

**Example**: Grid = (128, 128, 128), MPI = (3, 4, 2)

.. code-block:: text

   X dimension: 128 ÷ 3 = 42.666... ✗ NOT evenly divisible!
   Y dimension: 128 ÷ 4 = 32 ✓
   Z dimension: 128 ÷ 2 = 64 ✓
   
   Result: FAIL - Need different topology or grid size

Complete Implementation
------------------------

Here's the complete algorithm from HPC-ScaleTest:

.. code-block:: python

   # From engine/scaling.py
   
   class StrongScaling:
       """Generate strong scaling configurations.
       
       Strong scaling: Fixed problem size, increasing processes.
       Goal: Measure speedup as we add more computing resources.
       """
       
       def __init__(self, baseline_config, target_process_counts):
           """Initialize strong scaling generator.
           
           Args:
               baseline_config: Configuration with baseline process count
               target_process_counts: List of target process counts to test
           """
           self.baseline = baseline_config
           self.targets = target_process_counts
           
           # Extract baseline information
           self.baseline_procs = baseline_config.num_processes
           self.baseline_mpi = baseline_config.mpi_topology
           self.grid_dims = baseline_config.grid_dimensions
       
       def generate_configurations(self):
           """Generate all strong scaling configurations.
           
           Returns:
               List of configurations, one for each target process count
           """
           configs = []
           
           for target_procs in self.targets:
               # Calculate MPI topology for this process count
               mpi_topo = self._calculate_mpi_topology(target_procs)
               
               # Verify divisibility
               if not self._verify_divisibility(mpi_topo):
                   print(f"Warning: Process count {target_procs} incompatible with grid")
                   continue
               
               # Create configuration
               config = self._create_config(target_procs, mpi_topo)
               configs.append(config)
           
           return configs
       
       def _calculate_mpi_topology(self, num_procs):
           """Calculate balanced MPI topology.
           
           Uses prime factorization with balanced distribution.
           """
           if num_procs == 1:
               return (1, 1, 1)
           
           # Prime factorization
           factors = self._factorize(num_procs)
           
           # Distribute factors
           px, py, pz = 1, 1, 1
           for factor in sorted(factors, reverse=True):
               if px <= py and px <= pz:
                   px *= factor
               elif py <= pz:
                   py *= factor
               else:
                   pz *= factor
           
           return (px, py, pz)
       
       def _factorize(self, n):
           """Find all prime factors."""
           factors = []
           d = 2
           while d * d <= n:
               while n % d == 0:
                   factors.append(d)
                   n //= d
               d += 1
           if n > 1:
               factors.append(n)
           return factors
       
       def _verify_divisibility(self, mpi_topo):
           """Verify grid can be evenly divided."""
           nx, ny, nz = self.grid_dims
           px, py, pz = mpi_topo
           
           return (nx % px == 0 and 
                   ny % py == 0 and 
                   nz % pz == 0)
       
       def _create_config(self, num_procs, mpi_topo):
           """Create configuration for this process count."""
           config = {
               'num_processes': num_procs,
               'mpi_topology': mpi_topo,
               'grid_dimensions': self.grid_dims,  # UNCHANGED!
               'work_per_process': self._calculate_work_per_proc(num_procs)
           }
           return config
       
       def _calculate_work_per_proc(self, num_procs):
           """Calculate work per process."""
           nx, ny, nz = self.grid_dims
           total_work = nx * ny * nz
           return total_work / num_procs

Worked Example
===============

Complete Strong Scaling Study
-------------------------------

Let's walk through a complete example step-by-step.

**Scenario**: Climate simulation with 256 × 256 × 128 grid

Baseline Configuration
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   Application: Climate simulation
   Grid size: 256 × 256 × 128 cells (fixed for all runs!)
   Baseline: 8 processes on 1 node
   Target: Test 16, 32, 64, 128 processes

**Baseline MPI Topology**:

.. code-block:: text

   num_procs = 8
   factors = [2, 2, 2]  # 8 = 2³
   
   Distribution:
   - Factor 2: px = 2
   - Factor 2: py = 2
   - Factor 2: pz = 2
   
   Topology: (2, 2, 2)
   Verify: 2 × 2 × 2 = 8 ✓

**Baseline Work Distribution**:

.. code-block:: text

   Total cells = 256 × 256 × 128 = 8,388,608 cells
   
   Cells per process = 8,388,608 / 8 = 1,048,576 cells
   
   Per dimension:
   - X: 256 / 2 = 128 cells per process
   - Y: 256 / 2 = 128 cells per process
   - Z: 128 / 2 = 64 cells per process
   
   Each process handles: 128 × 128 × 64 = 1,048,576 cells ✓

Configuration 1: 16 Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Calculate MPI Topology**:

.. code-block:: text

   num_procs = 16
   factors = [2, 2, 2, 2]  # 16 = 2⁴
   
   Start: (1, 1, 1)
   
   Factor 2: px = 2 → (2, 1, 1)
   Factor 2: py = 2 → (2, 2, 1)
   Factor 2: pz = 2 → (2, 2, 2)
   Factor 2: px = 4 → (4, 2, 2)  # px was smallest
   
   Topology: (4, 2, 2)
   Verify: 4 × 2 × 2 = 16 ✓

**Verify Divisibility**:

.. code-block:: text

   Grid: (256, 256, 128)
   Topology: (4, 2, 2)
   
   X: 256 ÷ 4 = 64 ✓
   Y: 256 ÷ 2 = 128 ✓
   Z: 128 ÷ 2 = 64 ✓
   
   All evenly divisible! ✓

**Work Distribution**:

.. code-block:: text

   Total cells = 8,388,608 (unchanged)
   Processes = 16
   
   Cells per process = 8,388,608 / 16 = 524,288 cells
   
   Per dimension per process:
   - X: 256 / 4 = 64 cells
   - Y: 256 / 2 = 128 cells
   - Z: 128 / 2 = 64 cells
   
   Each process: 64 × 128 × 64 = 524,288 cells ✓
   
   Comparison to baseline:
   - Baseline: 1,048,576 cells per process
   - 16 procs: 524,288 cells per process
   - Reduction: 50% (half the work per process!)

**Expected Performance**:

.. code-block:: text

   If baseline takes 100 seconds:
   
   Ideal: 100 / 2 = 50 seconds (2× speedup from doubling processes)
   Realistic: ~55 seconds (communication overhead)
   
   Speedup: 100 / 55 = 1.82×
   Efficiency: 1.82 / 2 = 91%

Configuration 2: 32 Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**MPI Topology**:

.. code-block:: text

   num_procs = 32
   factors = [2, 2, 2, 2, 2]  # 32 = 2⁵
   
   Distribution → (4, 4, 2)
   Verify: 4 × 4 × 2 = 32 ✓

**Divisibility**:

.. code-block:: text

   256 ÷ 4 = 64 ✓
   256 ÷ 4 = 64 ✓
   128 ÷ 2 = 64 ✓

**Work Per Process**:

.. code-block:: text

   8,388,608 / 32 = 262,144 cells per process
   
   Per dimension: 64 × 64 × 64 = 262,144 ✓
   
   Reduction from baseline:
   1,048,576 → 262,144 (25% of original work)

Configuration 3: 64 Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**MPI Topology**:

.. code-block:: text

   num_procs = 64
   factors = [2, 2, 2, 2, 2, 2]  # 64 = 2⁶
   
   Distribution → (4, 4, 4)
   Verify: 4 × 4 × 4 = 64 ✓

**Work Per Process**:

.. code-block:: text

   8,388,608 / 64 = 131,072 cells per process
   
   Per dimension: 64 × 64 × 32 = 131,072 ✓
   
   Reduction: 12.5% of baseline work

Configuration 4: 128 Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**MPI Topology**:

.. code-block:: text

   num_procs = 128
   factors = [2, 2, 2, 2, 2, 2, 2]  # 128 = 2⁷
   
   Distribution → (8, 4, 4)
   Verify: 8 × 4 × 4 = 128 ✓

**Work Per Process**:

.. code-block:: text

   8,388,608 / 128 = 65,536 cells per process
   
   Per dimension: 32 × 64 × 32 = 65,536 ✓
   
   Reduction: 6.25% of baseline work

Summary Table
~~~~~~~~~~~~~~

.. code-block:: text

   ┌────────┬────────────────┬──────────────┬───────────────┬──────────┐
   │ Procs  │ MPI Topology   │ Work/Process │ Expected Time │ Speedup  │
   ├────────┼────────────────┼──────────────┼───────────────┼──────────┤
   │   8    │ (2, 2, 2)      │ 1,048,576    │ 100 sec       │ 1.0×     │
   │  16    │ (4, 2, 2)      │   524,288    │  55 sec       │ 1.82×    │
   │  32    │ (4, 4, 2)      │   262,144    │  30 sec       │ 3.33×    │
   │  64    │ (4, 4, 4)      │   131,072    │  17 sec       │ 5.88×    │
   │ 128    │ (8, 4, 4)      │    65,536    │  10 sec       │ 10.0×    │
   └────────┴────────────────┴──────────────┴───────────────┴──────────┘
   
   Note: Times are realistic estimates including communication overhead

Performance Analysis
=====================

Algorithmic Complexity
-----------------------

**Topology Calculation**: :math:`O(\log P)`

.. code-block:: text

   Prime factorization complexity:
   - Worst case: O(√P)
   - Typical case: O(log P) for powers of 2
   
   Factor distribution: O(log P)
   
   Total: O(log P) operations per configuration

**Example**:
- P = 1024: ~10 iterations
- P = 1,000,000: ~20 iterations
- Very fast even for large process counts!

**Complete Study**: :math:`O(N \log P_{\max})`

.. code-block:: text

   For N configurations testing up to Pmax processes:
   Time ≈ N × log(Pmax)
   
   Example: N = 10 configs, Pmax = 1024
   Time ≈ 10 × 10 = 100 operations
   
   On modern CPU: << 1 millisecond

Space Complexity
-----------------

**Per Configuration**: :math:`O(1)`

.. code-block:: text

   Each configuration stores:
   - Process count: 1 integer
   - MPI topology: 3 integers
   - Grid dimensions: 3 integers
   - Work metrics: 3 floats
   
   Total: ~50 bytes per configuration

**N Configurations**: :math:`O(N)`

.. code-block:: text

   10 configurations = 500 bytes
   100 configurations = 5 KB
   
   Negligible memory usage!

Communication Overhead
-----------------------

**Why Perfect Speedup is Impossible**

When processes communicate, they experience overhead:

1. **Latency**: Time to start sending a message (~1-10 microseconds)
2. **Bandwidth**: Time to transfer data (~1-10 GB/s)
3. **Synchronization**: Waiting for all processes to reach a point

**Communication Pattern in 3D Grid**:

.. code-block:: text

   Each process communicates with 6 neighbors (up/down/left/right/front/back)
   
   As we add more processes:
   - Work per process decreases (good!)
   - Communication per process increases (bad!)
   
   Eventually: Communication time > Computation time
   This limits maximum useful processes

**Example**: 256 × 256 × 128 grid

.. code-block:: text

   8 processes:
   - Work: 1,048,576 cells per process
   - Surface area: 3 × 128 × 64 = 24,576 cells to communicate
   - Communication/Computation: 2.3%
   
   128 processes:
   - Work: 65,536 cells per process
   - Surface area: 3 × 32 × 32 = 3,072 cells to communicate
   - Communication/Computation: 4.7%
   
   As processes increase, communication becomes more significant!

Efficiency Predictions
------------------------

**Theoretical Maximum**:

Based on Amdahl's Law with 95% parallelizable code:

.. code-block:: text

   Speedup(P) = P / (0.05 + 0.95/P)
   
   P=8:    Speedup = 7.27  (Efficiency = 91%)
   P=16:   Speedup = 13.91 (Efficiency = 87%)
   P=32:   Speedup = 24.62 (Efficiency = 77%)
   P=64:   Speedup = 39.02 (Efficiency = 61%)
   P=128:  Speedup = 58.18 (Efficiency = 45%)

**Key Insight**: Efficiency naturally decreases with more processes!

Practical Usage
================

Command-Line Interface
-----------------------

**Basic Strong Scaling Study**:

.. code-block:: bash

   # Generate strong scaling configurations
   hpc-scaletest strong-scaling \
       --baseline baseline_config.yaml \
       --nodes "2 4 8 16" \
       --procs-per-node 8
   
   # This creates configurations for:
   # - 16 processes (2 nodes × 8 procs)
   # - 32 processes (4 nodes × 8 procs)
   # - 64 processes (8 nodes × 8 procs)
   # - 128 processes (16 nodes × 8 procs)

**With Custom Process Counts**:

.. code-block:: bash

   # Test specific process counts
   hpc-scaletest strong-scaling \
       --baseline baseline_config.yaml \
       --process-counts "8 16 24 32 48 64"

**Output**:

.. code-block:: text

   Generated strong scaling configurations:
   ✓ config_8_procs.yaml   - MPI: (2,2,2)
   ✓ config_16_procs.yaml  - MPI: (4,2,2)
   ✓ config_24_procs.yaml  - MPI: (4,3,2)
   ✓ config_32_procs.yaml  - MPI: (4,4,2)
   ✓ config_48_procs.yaml  - MPI: (4,4,3)
   ✓ config_64_procs.yaml  - MPI: (4,4,4)
   
   All configurations verified for divisibility ✓
   
   Next steps:
   1. Review configurations in configs/
   2. Run: hpc-scaletest run --all
   3. Analyze: hpc-scaletest analyze results/

Python API Usage
-----------------

**Programmatic Generation**:

.. code-block:: text

   from hpc_scaletest import StrongScalingGenerator
   
   # Load baseline configuration
   baseline = {
       'num_processes': 8,
       'mpi_topology': (2, 2, 2),
       'grid_dimensions': (256, 256, 128)
   }
   
   # Create generator
   generator = StrongScalingGenerator(baseline)
   
   # Generate configurations for different process counts
   configs = generator.generate([16, 32, 64, 128])
   
   # Examine results
   for config in configs:
       print(f"Processes: {config.num_processes}")
       print(f"MPI Topology: {config.mpi_topology}")
       print(f"Work per process: {config.work_per_process:.0f} cells")
       print(f"Expected speedup: {config.expected_speedup:.2f}×")
       print()

**Output**:

.. code-block:: text

   Processes: 16
   MPI Topology: (4, 2, 2)
   Work per process: 524288 cells
   Expected speedup: 1.82×
   
   Processes: 32
   MPI Topology: (4, 4, 2)
   Work per process: 262144 cells
   Expected speedup: 3.33×
   
   Processes: 64
   MPI Topology: (4, 4, 4)
   Work per process: 131072 cells
   Expected speedup: 5.88×
   
   Processes: 128
   MPI Topology: (8, 4, 4)
   Work per process: 65536 cells
   Expected speedup: 10.00×

Advanced Usage
---------------

**Testing Different Grid Sizes**:

.. code-block:: text

   # Test if your grid size is suitable for strong scaling
   from hpc_scaletest import verify_strong_scaling
   
   grid = (256, 256, 128)
   process_counts = [8, 16, 32, 64, 128, 256]
   
   results = verify_strong_scaling(grid, process_counts)
   
   for procs, status in results.items():
       if status['divisible']:
           print(f"✓ {procs} processes: OK - MPI {status['topology']}")
       else:
           print(f"✗ {procs} processes: Cannot evenly divide grid")

**Custom MPI Topology Constraints**:

.. code-block:: python

   # Prefer certain dimensions (e.g., minimize Z-direction splitting)
   generator = StrongScalingGenerator(
       baseline,
       topology_preference='xy'  # Prefer splitting X and Y over Z
   )
   
   configs = generator.generate([64])
   # Result: Might get (8, 8, 1) instead of (4, 4, 4)

Testing and Validation
========================

Unit Tests
-----------

**Test 1: MPI Topology Calculation**:

.. code-block:: python

   def test_mpi_topology_calculation():
       """Verify MPI topology is correctly calculated."""
       
       # Test powers of 2
       assert calculate_mpi_topology(8) == (2, 2, 2)
       assert calculate_mpi_topology(64) == (4, 4, 4)
       
       # Test non-cubic numbers
       topo = calculate_mpi_topology(24)
       assert topo[0] * topo[1] * topo[2] == 24
       
       # Test prime numbers
       topo = calculate_mpi_topology(17)
       assert topo == (17, 1, 1) or topo == (1, 17, 1) or topo == (1, 1, 17)

**Test 2: Work Distribution**:

.. code-block:: python

   def test_work_distribution():
       """Verify work is correctly calculated per process."""
       
       baseline_work = 1_048_576  # cells per process at 8 procs
       grid = (256, 256, 128)
       
       # Test 16 processes
       work_16 = calculate_work_per_process(grid, 16)
       assert work_16 == 524_288
       assert work_16 == baseline_work / 2  # Half the work
       
       # Test 64 processes
       work_64 = calculate_work_per_process(grid, 64)
       assert work_64 == 131_072
       assert work_64 == baseline_work / 8  # 1/8 the work

**Test 3: Divisibility Verification**:

.. code-block:: python

   def test_divisibility():
       """Verify divisibility checks work correctly."""
       
       grid = (256, 256, 128)
       
       # Should pass
       assert verify_divisibility(grid, (4, 4, 2)) == True
       assert verify_divisibility(grid, (8, 8, 4)) == True
       
       # Should fail
       assert verify_divisibility(grid, (3, 3, 3)) == False  # 256 not divisible by 3
       assert verify_divisibility(grid, (5, 5, 5)) == False

Integration Tests
------------------

**End-to-End Test**:

.. code-block:: python

   def test_complete_strong_scaling_workflow():
       """Test complete workflow from input to output."""
       
       # Setup
       baseline = load_baseline_config("baseline.yaml")
       target_procs = [16, 32, 64]
       
       # Generate configurations
       generator = StrongScalingGenerator(baseline)
       configs = generator.generate(target_procs)
       
       # Verify all configs
       assert len(configs) == 3
       
       for config in configs:
           # Check MPI topology is valid
           px, py, pz = config.mpi_topology
           assert px * py * pz == config.num_processes
           
           # Check divisibility
           assert verify_divisibility(config.grid_dimensions, 
                                     config.mpi_topology)
           
           # Check work distribution
           expected_work = (baseline.total_cells / 
                           config.num_processes)
           assert abs(config.work_per_process - expected_work) < 1e-6

Verification Checklist
-----------------------

Before running your strong scaling study, verify:

.. code-block:: text

   ✓ Grid dimensions are fixed (same for all configs)
   ✓ Process counts increase (more processes = more speedup)
   ✓ All MPI topologies are valid (px × py × pz = total processes)
   ✓ Grid is evenly divisible for all process counts
   ✓ Work per process decreases proportionally
   ✓ No process has zero work
   ✓ Communication pattern is reasonable

**Automated Verification**:

.. code-block:: bash

   # Verify all configurations before running
   hpc-scaletest verify --type strong-scaling --configs configs/

Common Issues and Solutions
============================

Issue 1: Grid Not Divisible
-----------------------------

**Problem**:

.. code-block:: text

   Error: 256 cells not evenly divisible by 3 processes in X dimension

**Cause**: Your grid size (256) cannot be evenly divided by the number of processes (3).

**Solutions**:

1. **Change grid size** to be divisible by common process counts:

   .. code-block:: text
   
      Bad: 256 × 256 × 128 (not divisible by 3, 5, 6, 7, 9, ...)
      Good: 240 × 240 × 120 (divisible by 2, 3, 4, 5, 6, 8, ...)

2. **Use power-of-2 process counts**:

   .. code-block:: text
   
      Use: 8, 16, 32, 64, 128, 256, 512, 1024
      Avoid: 12, 18, 24, 36, 48 (unless grid matches)

3. **Adjust process count** to match grid:

   .. code-block:: bash
   
      # Find compatible process counts
      hpc-scaletest find-divisible-counts --grid 256 256 128
      # Output: 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096

Issue 2: Poor Speedup
-----------------------

**Problem**:

.. code-block:: text

   Expected: 8× faster with 64 processes
   Actual: Only 4× faster

**Causes**:

1. **Communication overhead**: Processes spend too much time communicating
2. **Load imbalance**: Some processes have more work than others
3. **Serial bottlenecks**: Parts of code that can't be parallelized
4. **I/O contention**: All processes writing to same disk

**Solutions**:

1. **Profile your application**:

   .. code-block:: bash
   
      # Use profiling tools to find bottlenecks
      mpirun -np 64 scalasca -analyze ./myapp
      # or
      mpirun -np 64 tau_exec ./myapp

2. **Optimize communication**:
   - Use non-blocking MPI calls
   - Reduce frequency of global synchronization
   - Overlap communication with computation

3. **Check load balance**:

   .. code-block:: python
   
      # Verify each process has equal work
      work_per_proc = [proc.cell_count for proc in all_processes]
      assert max(work_per_proc) - min(work_per_proc) < 100  # Small difference

Issue 3: Memory Limitations
-----------------------------

**Problem**:

.. code-block:: text

   Out of memory error with many processes

**Cause**: Each process needs its own memory copy of shared data.

**Solutions**:

1. **Reduce per-process memory**:
   - More processes = less data per process
   - Strong scaling naturally helps with memory!

2. **Check memory usage**:

   .. code-block:: bash
   
      # Monitor memory per process
      mpirun -np 64 /usr/bin/time -v ./myapp
      # Look for "Maximum resident set size"

3. **Use memory-efficient algorithms**:
   - Streaming algorithms
   - Out-of-core methods
   - Data compression

Best Practices
===============

Choosing Process Counts
-------------------------

**Rule 1: Start with powers of 2**

.. code-block:: text

   Good: 1, 2, 4, 8, 16, 32, 64, 128, 256
   
   Why: Most HPC systems have power-of-2 cores per node
        Most grid sizes are powers of 2
        Simplifies MPI topology

**Rule 2: Match your node architecture**

.. code-block:: text

   Example: 8 cores per node
   
   Test: 1, 2, 4 nodes (8, 16, 32 processes)
   
   Example: 36 cores per node
   
   Test: 1, 2, 4, 8 nodes (36, 72, 144, 288 processes)

**Rule 3: Stop when efficiency drops below 50%**

.. code-block:: text

   If you're getting less than 50% efficiency:
   - Communication overhead is too high
   - Problem is too small for that many processes
   - Time to consider weak scaling instead

Designing Grid Sizes
----------------------

**Highly divisible numbers are best**:

.. code-block:: text

   Excellent: 128, 256, 512, 1024 (powers of 2)
   Good: 120, 240, 360, 480 (divisible by many numbers)
   Poor: 127, 251, 509 (prime numbers)

**Calculator for grid sizes**:

.. code-block:: python

   def suggest_grid_size(approx_size, max_procs):
       """Suggest grid size divisible by common process counts."""
       
       # Find highly composite number near approx_size
       factors_needed = factorize(max_procs)
       
       # Adjust size to be divisible by all factors
       adjusted = approx_size
       for factor in set(factors_needed):
           while adjusted % factor != 0:
               adjusted += 1
       
       return adjusted
   
   # Example
   size = suggest_grid_size(250, 64)
   print(size)  # 256 (divisible by 1,2,4,8,16,32,64)

Running Studies Efficiently
-----------------------------

**1. Test locally first**:

.. code-block:: bash

   # Run small test on your laptop
   hpc-scaletest strong-scaling \
       --baseline tiny_baseline.yaml \
       --process-counts "1 2 4" \
       --dry-run

**2. Start with small process counts**:

.. code-block:: bash

   # First study: just 2× and 4× scaling
   hpc-scaletest strong-scaling \
       --baseline baseline.yaml \
       --process-counts "8 16 32"

**3. Scale up gradually**:

.. code-block:: bash

   # If 32 processes worked well, try more
   hpc-scaletest strong-scaling \
       --baseline baseline.yaml \
       --process-counts "32 64 128 256"

**4. Analyze after each batch**:

.. code-block:: bash

   # Check results before running more
   hpc-scaletest analyze results/
   # Look at efficiency plot
   # If efficiency < 50%, stop adding processes

Summary
========

Key Takeaways
--------------

**What is Strong Scaling?**

- Fixed problem size
- Increasing number of processes
- Goal: Faster results

**How Does HPC-ScaleTest Help?**

- Automatically calculates MPI topologies
- Verifies grid divisibility
- Generates all configuration files
- Ensures proper work distribution

**Important Concepts**

1. **Speedup**: How much faster with N processes
2. **Efficiency**: How well you're using the processes
3. **Amdahl's Law**: Serial parts limit maximum speedup
4. **Communication Overhead**: Limits practical process count

**Best Practices**

- Use power-of-2 process counts
- Choose highly divisible grid sizes
- Stop when efficiency drops below 50%
- Profile to find bottlenecks

Further Reading
----------------

**Parallel Computing Concepts**:
- "Introduction to Parallel Computing" by Grama et al.
- "Parallel Programming in C with MPI and OpenMP" by Quinn

**Performance Analysis**:
- "Performance Optimization of Numerically Intensive Codes" by Goedecker & Hoisie
- "Tools and Technologies for Performance Analysis" by various authors

**HPC-ScaleTest Documentation**:
- :doc:`weak_scaling_complete` - Weak scaling algorithm
- :doc:`mpi_topology_complete` - MPI topology details
- :doc:`../user_guide` - Complete user guide

Mathematical Appendix
======================

Speedup Formula Derivation
----------------------------

**Given**:

.. math::

   T_1 = \text{Time on 1 process}

.. math::

   T_P = \text{Time on } P \text{ processes}

**Speedup Definition**:

.. math::

   S(P) = \frac{T_1}{T_P}

**Amdahl's Law**:

Let :math:`f` be the fraction of serial code (0 ≤ f ≤ 1):

.. math::

   T_1 = f \cdot T_1 + (1-f) \cdot T_1

.. math::

   T_P = f \cdot T_1 + \frac{(1-f) \cdot T_1}{P}

**Speedup**:

.. math::

   S(P) = \frac{T_1}{f \cdot T_1 + \frac{(1-f) \cdot T_1}{P}}

.. math::

   S(P) = \frac{1}{f + \frac{1-f}{P}}

**Maximum Speedup** (as :math:`P \to \infty`):

.. math::

   S_{\max} = \lim_{P \to \infty} S(P) = \frac{1}{f}

**Example**: If :math:`f = 0.05` (5% serial), then :math:`S_{\max} = 20`.

Efficiency Formula
-------------------

**Definition**:

.. math::

   E(P) = \frac{S(P)}{P} = \frac{T_1}{P \cdot T_P}

**Perfect Efficiency**: :math:`E(P) = 1` (100%)

**Typical**: :math:`E(P)` decreases as :math:`P` increases

**From Amdahl's Law**:

.. math::

   E(P) = \frac{1}{P \left( f + \frac{1-f}{P} \right)} = \frac{1}{1 + f(P-1)}

Work Distribution Formula
---------------------------

**Total Work**:

.. math::

   W = n_x \times n_y \times n_z

**Work Per Process**:

.. math::

   W_P = \frac{W}{P} = \frac{n_x \times n_y \times n_z}{p_x \times p_y \times p_z}

**Per Dimension**:

.. math::

   w_x = \frac{n_x}{p_x}, \quad w_y = \frac{n_y}{p_y}, \quad w_z = \frac{n_z}{p_z}

**Verification**:

.. math::

   W_P = w_x \times w_y \times w_z = \frac{n_x}{p_x} \times \frac{n_y}{p_y} \times \frac{n_z}{p_z} = \frac{W}{P} \quad \checkmark
