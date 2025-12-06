================================================================
Weak Scaling Algorithm: Complete Technical Reference
================================================================

.. contents:: Table of Contents
   :local:
   :depth: 3

Introduction for Non-Technical Readers
========================================

What is This Document About?
-----------------------------

This document explains **how HPC-ScaleTest automatically calculates the right settings** 
for running your scientific application on multiple computers (called "nodes") at the same time.

Think of it like this:

**Scenario**: You have a weather simulation program that takes 1 hour to run on 1 computer.
You want to simulate a bigger area (2x bigger), but you don't want it to take 2 hours.

**Solution**: Use 2 computers at the same time! HPC-ScaleTest figures out exactly how to 
split the work so each computer still takes about 1 hour.

**That's "weak scaling"** - making the problem bigger while keeping the time the same by 
using more computers.

Who Should Read This?
----------------------

This document is for:

1. **Users** - Learn how HPC-ScaleTest works (no coding required!)
2. **Developers** - Understand the algorithms to modify or extend them
3. **Researchers** - See the mathematical foundations for publications
4. **Students** - Study parallel computing with real-world examples

You don't need programming experience to understand the concepts - we explain everything!

What You'll Learn
-----------------

After reading this document, you'll understand:

✅ **What weak scaling means** (in plain English)
✅ **Why it's useful** (save time on big simulations!)
✅ **How HPC-ScaleTest calculates it** (the math, explained simply)
✅ **How to use it** (practical examples)
✅ **How to verify it works** (check your results)

Reading Guide
-------------

**If you're new to HPC**:
   Start with "Conceptual Foundation" - it explains everything from scratch.

**If you're a developer**:
   Jump to "Algorithm Design" for implementation details.

**If you're a researcher**:
   See "Mathematical Foundations" for formal definitions and proofs.

**If you just want to use it**:
   Skip to "Practical Usage" for command-line examples.

**A comprehensive deep-dive into the weak scaling algorithm implementation**

.. contents:: Table of Contents
   :local:
   :depth: 4

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PART 1: MATHEMATICAL FOUNDATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1.1 Weak Scaling Definition
============================

Weak scaling maintains **constant work per processor** while increasing the total problem size proportionally.

Mathematical Formulation
------------------------

Given:
- W = Work (computational load)
- P = Number of processors
- T = Execution time

**Strong Scaling**: W = constant, P increases → T should decrease proportionally
**Weak Scaling**: W/P = constant, both W and P increase → T should stay constant

For grid-based simulations:

.. math::

   W \propto N_{cells} \times N_{timesteps}

In weak scaling:

.. math::

   \frac{N_{cells}}{P} = \text{constant}

This means:

.. math::

   N_{cells}(P) = N_{cells}(P_0) \times \frac{P}{P_0}

Where P₀ is the baseline processor count.

1.2 Domain Decomposition Mathematics
====================================

3D Grid Decomposition
---------------------

A 3D computational domain is decomposed across P processors in three dimensions:

.. math::

   P = P_x \times P_y \times P_z

Where Pₓ, Pᵧ, Pᵤ are the number of processes in each dimension.

The global grid size is:

.. math::

   N_{total} = N_x \times N_y \times N_z

Each processor owns a subdomain:

.. math::

   n_{local} = \frac{N_x}{P_x} \times \frac{N_y}{P_y} \times \frac{N_z}{P_z}

For weak scaling, n_local must remain constant:

.. math::

   n_{local} = \frac{N_x(P)}{P_x(P)} \times \frac{N_y(P)}{P_y(P)} \times \frac{N_z(P)}{P_z(P)} = \text{constant}

Physical Domain Scaling
------------------------

The physical domain size L also scales proportionally:

.. math::

   \begin{aligned}
   L_x(P) &= L_x(P_0) \times \frac{P_x(P)}{P_x(P_0)} \\
   L_y(P) &= L_y(P_0) \times \frac{P_y(P)}{P_y(P_0)} \\
   L_z(P) &= L_z(P_0) \times \frac{P_z(P)}{P_z(P_0)}
   \end{aligned}

This ensures grid resolution (Δx, Δy, Δz) remains constant:

.. math::

   \Delta x = \frac{L_x}{N_x} = \text{constant}

1.3 Scaling Factor Mathematics
==============================

Definition
----------

A scaling factor s determines how aggressively the problem grows:

.. math::

   P_i = P_0 \times s^i

where i is the scaling step number.

Common choices:
- s = 2 (doubling, most common)
- s = 4 (quadrupling)
- s = √2 (fine-grained)

For doubling (s = 2):

.. code-block:: text

   Step 0: P₀ processors
   Step 1: 2P₀ processors (doubled)
   Step 2: 4P₀ processors (quadrupled)
   Step 3: 8P₀ processors
   Step i: 2ⁱP₀ processors

Dimension-Specific Scaling
---------------------------

When scaling in dimension d with factor s:

.. math::

   P_d^{new} = P_d^{old} \times s

The corresponding grid and domain also scale:

.. math::

   \begin{aligned}
   N_d^{new} &= N_d^{old} \times s \\
   L_d^{new} &= L_d^{old} \times s
   \end{aligned}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PART 2: ALGORITHM DESIGN
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

2.1 Dimension Rotation Strategy
===============================

The Problem
-----------

Naïve weak scaling might scale all dimensions equally:

.. code-block:: text

   # WRONG: Isotropic scaling
   Pₓ, Pᵧ, Pᵤ = 2, 2, 2  # Start
   # Scale all equally
   Pₓ, Pᵧ, Pᵤ = 4, 4, 4  # 64× increase (too aggressive!)

This causes **exponential growth** and is impractical.

The Solution: Rotation
----------------------

Scale ONE dimension at a time in a rotating pattern:

**3D Pattern** (X → Y → Z → X → Y → Z ...):

.. code-block:: text

   i=0: (2,2,2) → 8 procs   [baseline]
   i=1: (4,2,2) → 16 procs  [scale X]
   i=2: (4,4,2) → 32 procs  [scale Y]
   i=3: (4,4,4) → 64 procs  [scale Z]
   i=4: (8,4,4) → 128 procs [scale X again]
   i=5: (8,8,4) → 256 procs [scale Y again]

**2D Pattern** (X → Y → X → Y ...):

.. code-block:: text

   i=0: (2,2,2) → 8 procs   [baseline, Z fixed]
   i=1: (4,2,2) → 16 procs  [scale X]
   i=2: (4,4,2) → 32 procs  [scale Y]
   i=3: (8,4,2) → 64 procs  [scale X]
   i=4: (8,8,2) → 128 procs [scale Y]

**1D Pattern** (X only):

.. code-block:: text

   i=0: (2,2,2) → 8 procs   [baseline, Y,Z fixed]
   i=1: (4,2,2) → 16 procs  [scale X]
   i=2: (8,2,2) → 32 procs  [scale X]
   i=3: (16,2,2) → 64 procs [scale X]

Mathematical Formulation
------------------------

For 3D scaling with factor s=2:

.. math::

   \begin{aligned}
   P_x(i) &= P_x(0) \times 2^{\lfloor (i+2)/3 \rfloor} \\
   P_y(i) &= P_y(0) \times 2^{\lfloor (i+1)/3 \rfloor} \\
   P_z(i) &= P_z(0) \times 2^{\lfloor i/3 \rfloor}
   \end{aligned}

More generally, for dimension d at step i:

.. math::

   P_d(i) = P_d(0) \times s^{n_d(i)}

where n_d(i) is the number of times dimension d has been scaled up to step i.

For 3D with rotation:

.. code-block:: python

   def num_scalings_for_dim(step, dim, total_dims=3):
       """How many times has this dimension been scaled?"""
       # dim: 0=X, 1=Y, 2=Z
       # Rotation pattern: X(0), Y(1), Z(2), X(0), Y(1), Z(2), ...
       
       full_cycles = step // total_dims
       remainder = step % total_dims
       
       if dim < remainder:
           return full_cycles + 1
       else:
           return full_cycles

2.2 Implementation Algorithm
============================

Pseudocode
----------

.. code-block:: python

   def compute_weak_scaling_configuration(
       baseline,        # Baseline configuration (node 1)
       target_procs,    # Target number of processes
       scale_factor,    # Scaling factor (e.g., 2)
       dimensions       # 1, 2, or 3
   ):
       """
       Compute weak scaling configuration for target process count.
       
       Args:
           baseline: dict with 'Px', 'Py', 'Pz', 'Lx', 'Ly', 'Lz', 'Nx', 'Ny', 'Nz'
           target_procs: Target total processes
           scale_factor: Scaling factor per step
           dimensions: Number of dimensions to scale (1, 2, or 3)
       
       Returns:
           dict: Scaled configuration
       """
       
       # Step 1: Calculate number of scaling steps
       baseline_procs = baseline['Px'] * baseline['Py'] * baseline['Pz']
       ratio = target_procs / baseline_procs
       steps = log_base(ratio, scale_factor)
       
       # Step 2: Initialize from baseline
       config = baseline.copy()
       
       # Step 3: Apply rotational scaling
       for step in range(steps):
           # Determine which dimension to scale
           dim = step % dimensions  # 0=X, 1=Y, 2=Z
           
           if dim == 0:  # Scale X
               config['Px'] *= scale_factor
               config['Lx'] *= scale_factor
               config['Nx'] *= scale_factor
           elif dim == 1:  # Scale Y
               config['Py'] *= scale_factor
               config['Ly'] *= scale_factor
               config['Ny'] *= scale_factor
           elif dim == 2:  # Scale Z (only if 3D)
               config['Pz'] *= scale_factor
               config['Lz'] *= scale_factor
               config['Nz'] *= scale_factor
       
       return config

Detailed Implementation with Real Code
---------------------------------------

From ``engine/scaling.py``:

.. code-block:: text
   :linenos:
   :emphasize-lines: 15-17,29-31

   def _compute_mpi_topology_3d(self, node_count, target_procs):
       """
       Compute 3D MPI topology using dimension rotation.
       
       This implements the core weak scaling algorithm for 3D problems.
       Dimensions are scaled in rotation (X→Y→Z→X...) to maintain
       balanced decomposition while achieving target process count.
       
       Args:
           node_count: Current node number (for logging)
           target_procs: Target total number of MPI processes
       
       Returns:
           tuple: (nproc_x, nproc_y, nproc_z) MPI topology
       """
       # Extract baseline configuration
       baseline_x = self.base_config['nproc_x']
       baseline_y = self.base_config['nproc_y']
       baseline_z = self.base_config['nproc_z']
       baseline_procs = baseline_x * baseline_y * baseline_z
       
       # Calculate scaling ratio and steps
       if target_procs == baseline_procs:
           # Node 1: return baseline exactly
           return (baseline_x, baseline_y, baseline_z)
       
       # Calculate number of doublings (assuming scale_factor=2)
       ratio = target_procs / baseline_procs
       doublings = int(round(math.log2(ratio)))
       
       logger.debug(f"  3D Scaling: {baseline_procs} → {target_procs}")
       logger.debug(f"    Ratio: {ratio}, Doublings: {doublings}")
       
       # Apply rotational scaling
       nproc_x, nproc_y, nproc_z = baseline_x, baseline_y, baseline_z
       
       for i in range(doublings):
           # Determine dimension to scale: i % 3
           # 0 → X, 1 → Y, 2 → Z, 3 → X, ...
           dim_index = i % 3
           
           logger.debug(f"    Step {i}: dim={dim_index} ", end="")
           
           if dim_index == 0:  # Scale X
               nproc_x *= self.scale_factor
               logger.debug(f"X: {nproc_x//self.scale_factor}→{nproc_x}")
           elif dim_index == 1:  # Scale Y
               nproc_y *= self.scale_factor
               logger.debug(f"Y: {nproc_y//self.scale_factor}→{nproc_y}")
           elif dim_index == 2:  # Scale Z
               nproc_z *= self.scale_factor
               logger.debug(f"Z: {nproc_z//self.scale_factor}→{nproc_z}")
       
       # Verify we hit the target
       final_procs = nproc_x * nproc_y * nproc_z
       if final_procs != target_procs:
           logger.warning(
               f"  MPI topology mismatch: "
               f"computed {final_procs} != target {target_procs}"
           )
       
       return (int(nproc_x), int(nproc_y), int(nproc_z))

2.3 Data Structure Design
=========================

Configuration Dictionary
------------------------

The algorithm uses a configuration dictionary:

.. code-block:: python

   config = {
       # Node information
       'nodes': 4,              # Number of compute nodes
       'total_procs': 32,       # Total MPI processes
       
       # MPI topology
       'nproc_x': 4,            # Processes in X dimension
       'nproc_y': 4,            # Processes in Y dimension
       'nproc_z': 2,            # Processes in Z dimension
       
       # Physical domain (simulation space)
       'Lx': 20.0,              # Domain length in X (meters, km, etc.)
       'Ly': 20.0,              # Domain length in Y
       'Lz': 10.0,              # Domain length in Z
       
       # Computational grid
       'nx': 256,               # Grid cells in X
       'ny': 256,               # Grid cells in Y
       'nz': 128,               # Grid cells in Z
       
       # Particles per cell (for particle codes)
       'num_particles_x': 100,  # Particles in X per cell
       'num_particles_y': 100,  # Particles in Y per cell
       'num_particles_z': 100   # Particles in Z per cell
   }

Invariants
----------

The configuration must satisfy these invariants:

1. **Process Count Match**:

   .. code-block:: python
   
      assert config['total_procs'] == (
          config['nproc_x'] * 
          config['nproc_y'] * 
          config['nproc_z']
      )

2. **Constant Work Per Process**:

   .. code-block:: python
   
      cells_per_proc = (
          config['nx'] * config['ny'] * config['nz']
      ) / config['total_procs']
      
      # Must be constant across all configurations
      assert abs(cells_per_proc - BASELINE_CELLS_PER_PROC) < 0.01

3. **Constant Grid Resolution**:

   .. code-block:: python
   
      delta_x = config['Lx'] / config['nx']
      delta_y = config['Ly'] / config['ny']
      delta_z = config['Lz'] / config['nz']
      
      # Must match baseline resolution
      assert abs(delta_x - BASELINE_DELTA_X) < 1e-10

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PART 3: WORKED EXAMPLE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

3.1 Complete 3D Weak Scaling Example
====================================

Problem Setup
-------------

**Hardware Configuration**:

.. code-block:: yaml

   hardware:
     procs_per_node: 8    # 8 MPI ranks per node
     nodes: [1, 2, 4, 8]  # Test these node counts

**Initial Configuration** (Node 1):

.. code-block:: yaml

   scaling:
     initial_procs: [2, 2, 2]          # 2×2×2 = 8 processes
     initial_domain: [10.0, 10.0, 10.0] # 10×10×10 physical domain
     initial_cells: [128, 128, 128]     # 128³ computational cells
     scale_factor: 2                    # Double each step
     scaling_dimensions: 3              # 3D scaling

Node 1: Baseline Configuration
-------------------------------

**Input**:

.. code-block:: text

   node_count = 1
   target_procs = 1 × 8 = 8

**Calculation**:

.. code-block:: python

   # This is the baseline - no scaling applied
   config_n1 = {
       'nodes': 1,
       'total_procs': 8,
       
       # MPI topology (from YAML)
       'nproc_x': 2,
       'nproc_y': 2,
       'nproc_z': 2,
       
       # Domain (from YAML)
       'Lx': 10.0,
       'Ly': 10.0,
       'Lz': 10.0,
       
       # Grid (from YAML)
       'nx': 128,
       'ny': 128,
       'nz': 128,
   }

**Verification**:

.. code-block:: text

   # Total processes
   total = 2 × 2 × 2 = 8 ✓
   
   # Cells per process
   total_cells = 128 × 128 × 128 = 2,097,152
   cells_per_proc = 2,097,152 / 8 = 262,144 cells/process
   
   # Cells per process per dimension
   nx_per_proc = 128 / 2 = 64
   ny_per_proc = 128 / 2 = 64
   nz_per_proc = 128 / 2 = 64
   
   # Grid resolution
   Δx = 10.0 / 128 = 0.078125
   Δy = 10.0 / 128 = 0.078125
   Δz = 10.0 / 128 = 0.078125

Node 2: First Scaling Step (X)
-------------------------------

**Input**:

.. code-block:: text

   node_count = 2
   target_procs = 2 × 8 = 16

**Step-by-Step Calculation**:

.. code-block:: text

   # Step 1: Calculate scaling
   baseline_procs = 8
   ratio = 16 / 8 = 2
   doublings = log2(2) = 1  # One doubling step
   
   # Step 2: Apply rotational scaling
   # Starting point
   nproc_x, nproc_y, nproc_z = 2, 2, 2
   
   # Loop: i = 0 (first doubling)
   i = 0
   dim_index = 0 % 3 = 0  # Scale X dimension
   
   # Scale X
   nproc_x *= 2  # 2 → 4
   # nproc_y unchanged (still 2)
   # nproc_z unchanged (still 2)
   
   # Result after loop
   nproc_x, nproc_y, nproc_z = 4, 2, 2
   
   # Step 3: Calculate scaling factors
   scale_x = 4 / 2 = 2.0  # X doubled
   scale_y = 2 / 2 = 1.0  # Y unchanged
   scale_z = 2 / 2 = 1.0  # Z unchanged
   
   # Step 4: Scale domain
   Lx = 10.0 × 2.0 = 20.0  # Doubled
   Ly = 10.0 × 1.0 = 10.0  # Unchanged
   Lz = 10.0 × 1.0 = 10.0  # Unchanged
   
   # Step 5: Scale grid
   nx = 128 × 2 = 256  # Doubled
   ny = 128 × 1 = 128  # Unchanged
   nz = 128 × 1 = 128  # Unchanged

**Result**:

.. code-block:: python

   config_n2 = {
       'nodes': 2,
       'total_procs': 16,
       
       'nproc_x': 4,   # ← Scaled
       'nproc_y': 2,   # ← Fixed
       'nproc_z': 2,   # ← Fixed
       
       'Lx': 20.0,     # ← Scaled
       'Ly': 10.0,     # ← Fixed
       'Lz': 10.0,     # ← Fixed
       
       'nx': 256,      # ← Scaled
       'ny': 128,      # ← Fixed
       'nz': 128,      # ← Fixed
   }

**Verification**:

.. code-block:: text

   # Total processes
   total = 4 × 2 × 2 = 16 ✓
   
   # Cells per process
   total_cells = 256 × 128 × 128 = 4,194,304
   cells_per_proc = 4,194,304 / 16 = 262,144 ✓
   # Same as Node 1!
   
   # Cells per process per dimension
   nx_per_proc = 256 / 4 = 64 ✓
   ny_per_proc = 128 / 2 = 64 ✓
   nz_per_proc = 128 / 2 = 64 ✓
   # Same distribution as Node 1!
   
   # Grid resolution
   Δx = 20.0 / 256 = 0.078125 ✓
   Δy = 10.0 / 128 = 0.078125 ✓
   Δz = 10.0 / 128 = 0.078125 ✓
   # Same resolution as Node 1!

Node 4: Second Scaling Step (Y)
--------------------------------

**Input**:

.. code-block:: text

   node_count = 4
   target_procs = 4 × 8 = 32

**Step-by-Step Calculation**:

.. code-block:: text

   # Step 1: Calculate scaling
   baseline_procs = 8
   ratio = 32 / 8 = 4
   doublings = log2(4) = 2  # Two doubling steps
   
   # Step 2: Apply rotational scaling
   nproc_x, nproc_y, nproc_z = 2, 2, 2  # Start from baseline
   
   # Loop iteration 1: i = 0
   dim_index = 0 % 3 = 0  # Scale X
   nproc_x *= 2  # 2 → 4
   # After i=0: (4, 2, 2)
   
   # Loop iteration 2: i = 1
   dim_index = 1 % 3 = 1  # Scale Y
   nproc_y *= 2  # 2 → 4
   # After i=1: (4, 4, 2)
   
   # Result
   nproc_x, nproc_y, nproc_z = 4, 4, 2
   
   # Step 3: Scaling factors
   scale_x = 4 / 2 = 2.0
   scale_y = 4 / 2 = 2.0  # Now Y also scales
   scale_z = 2 / 2 = 1.0
   
   # Step 4: Scale domain
   Lx = 10.0 × 2.0 = 20.0
   Ly = 10.0 × 2.0 = 20.0  # Now Y scales
   Lz = 10.0 × 1.0 = 10.0
   
   # Step 5: Scale grid
   nx = 128 × 2 = 256
   ny = 128 × 2 = 256  # Now Y scales
   nz = 128 × 1 = 128

**Result**:

.. code-block:: python

   config_n4 = {
       'nodes': 4,
       'total_procs': 32,
       
       'nproc_x': 4,   # ← Scaled (step 0)
       'nproc_y': 4,   # ← Scaled (step 1)
       'nproc_z': 2,   # ← Still fixed
       
       'Lx': 20.0,     # ← Scaled
       'Ly': 20.0,     # ← Scaled
       'Lz': 10.0,     # ← Still fixed
       
       'nx': 256,      # ← Scaled
       'ny': 256,      # ← Scaled
       'nz': 128,      # ← Still fixed
   }

**Verification**:

.. code-block:: text

   # Total processes
   total = 4 × 4 × 2 = 32 ✓
   
   # Cells per process
   total_cells = 256 × 256 × 128 = 8,388,608
   cells_per_proc = 8,388,608 / 32 = 262,144 ✓
   # STILL same as baseline!
   
   # Per-dimension distribution
   nx_per_proc = 256 / 4 = 64 ✓
   ny_per_proc = 256 / 4 = 64 ✓
   nz_per_proc = 128 / 2 = 64 ✓
   # STILL (64, 64, 64) per process!
   
   # Grid resolution
   Δx = 20.0 / 256 = 0.078125 ✓
   Δy = 20.0 / 256 = 0.078125 ✓
   Δz = 10.0 / 128 = 0.078125 ✓
   # STILL constant!

[Content continues with Node 8 calculation and complete verification...]

3.2 Complete Configuration Table
================================

.. code-block:: text

   ╔══════╦═════════╦═══════════════════╦════════════════════════╦═══════════════════════╗
   ║ Node ║  Procs  ║  MPI Topology     ║  Domain Size           ║  Grid Cells           ║
   ║  #   ║  Total  ║  (Px × Py × Pz)   ║  (Lx, Ly, Lz)          ║  (nx, ny, nz)         ║
   ╠══════╬═════════╬═══════════════════╬════════════════════════╬═══════════════════════╣
   ║  1   ║    8    ║   2 ×  2 ×  2     ║  10.0, 10.0, 10.0      ║  128,  128,  128      ║
   ║      ║         ║                   ║  Δ=0.078125            ║  64×64×64 per proc    ║
   ╠══════╬═════════╬═══════════════════╬════════════════════════╬═══════════════════════╣
   ║  2   ║   16    ║   4 ×  2 ×  2     ║  20.0, 10.0, 10.0      ║  256,  128,  128      ║
   ║      ║         ║   [X scaled]      ║  [X doubled]           ║  [X doubled]          ║
   ║      ║         ║                   ║  Δ=0.078125 ✓          ║  64×64×64 per proc ✓  ║
   ╠══════╬═════════╬═══════════════════╬════════════════════════╬═══════════════════════╣
   ║  4   ║   32    ║   4 ×  4 ×  2     ║  20.0, 20.0, 10.0      ║  256,  256,  128      ║
   ║      ║         ║   [X,Y scaled]    ║  [X,Y doubled]         ║  [X,Y doubled]        ║
   ║      ║         ║                   ║  Δ=0.078125 ✓          ║  64×64×64 per proc ✓  ║
   ╠══════╬═════════╬═══════════════════╬════════════════════════╬═══════════════════════╣
   ║  8   ║   64    ║   4 ×  4 ×  4     ║  20.0, 20.0, 20.0      ║  256,  256,  256      ║
   ║      ║         ║   [X,Y,Z scaled]  ║  [X,Y,Z doubled]       ║  [X,Y,Z doubled]      ║
   ║      ║         ║                   ║  Δ=0.078125 ✓          ║  64×64×64 per proc ✓  ║
   ╠══════╬═════════╬═══════════════════╬════════════════════════╬═══════════════════════╣
   ║  16  ║  128    ║   8 ×  4 ×  4     ║  40.0, 20.0, 20.0      ║  512,  256,  256      ║
   ║      ║         ║   [X scaled 2×]   ║  [X quadrupled]        ║  [X quadrupled]       ║
   ║      ║         ║                   ║  Δ=0.078125 ✓          ║  64×64×64 per proc ✓  ║
   ╚══════╩═════════╩═══════════════════╩════════════════════════╩═══════════════════════╝

**Observations**:

1. Grid resolution Δ **constant** ✓
2. Cells per process **constant** (64×64×64 = 262,144) ✓
3. MPI topology follows rotation pattern (X→Y→Z→X) ✓
4. Domain scales proportionally with processes ✓

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PART 4: PERFORMANCE ANALYSIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

4.1 Algorithmic Complexity
==========================

Time Complexity
---------------

For computing a single configuration:

**MPI Topology Calculation**:

.. code-block:: text

   # Number of iterations = log₂(target_procs / baseline_procs)
   T_topology = O(log₂(P/P₀))
   
   # For P = 512, P₀ = 8:
   T_topology = O(log₂(64)) = O(6) iterations

**Domain Scaling**:

.. code-block:: python

   # Three multiplications (Lx, Ly, Lz)
   T_domain = O(1)

**Grid Scaling**:

.. code-block:: python

   # Three multiplications (nx, ny, nz)
   T_grid = O(1)

**Input File Generation**:

.. code-block:: text

   # Regex replacements = O(M × K)
   # M = input file size, K = number of parameters
   T_input = O(M × K)
   
   # Typically: M ~ 1000 lines, K ~ 10 parameters
   T_input = O(10,000) operations

**Total Per Configuration**:

.. code-block:: text

   T_total = O(log P + M×K)
   
   # Usually dominated by input file generation
   T_total ≈ O(M×K)

For N node configurations:

.. code-block:: text

   T_complete = N × O(log P + M×K)
   
   # Example: N=7 configurations, M=1000, K=10
   T_complete ≈ 7 × 10,000 = 70,000 operations
   
   # On modern CPU: < 1 second

Space Complexity
----------------

**Configuration Storage**:

.. code-block:: text

   # Each config: ~15 numeric values
   S_config = O(1) ≈ 15 × 8 bytes = 120 bytes
   
   # N configurations
   S_total = N × 120 bytes
   
   # Example: 10 configs = 1.2 KB (negligible)

**Input File Storage**:

.. code-block:: text

   # M lines × ~80 chars/line
   S_input = O(M) ≈ M × 80 bytes
   
   # Example: 1000 lines = 80 KB per file
   # N files = N × 80 KB
   # 10 files = 800 KB (negligible)

**Conclusion**: Memory usage is **negligible** for typical use cases.

4.2 Numerical Stability
=======================

Floating Point Considerations
------------------------------

The algorithm uses floating point arithmetic for domain sizes:

.. code-block:: text

   Lx_new = Lx_baseline × scale_factor^i

**Potential Issue**: Accumulated rounding errors

**Example**:

.. code-block:: text

   # Baseline
   Lx = 10.0
   scale_factor = 2.0
   
   # After 10 doublings
   Lx_10 = 10.0 × (2.0 ** 10) = 10240.0
   
   # Floating point representation
   actual = 10240.000000000000  # Exact in float64
   
   # Error: 0.0 ✓ (no error for powers of 2)

**However**, for non-power-of-2 factors:

.. code-block:: text

   # scale_factor = 1.5
   Lx_10 = 10.0 × (1.5 ** 10) = 576.650390625
   
   # Actual float64
   actual ≈ 576.6503906250001
   
   # Relative error
   error = |576.650390625 - 576.6503906250001| / 576.650390625
       ≈ 1.7e-16 (negligible)

Integer Arithmetic
------------------

Grid cell counts use integer arithmetic:

.. code-block:: text

   nx_new = nx_baseline × scale_factor^i

**Guaranteed Exact** for integer scale factors (2, 4, etc.):

.. code-block:: text

   nx = 128
   scale = 2
   
   # After 10 doublings
   nx_10 = 128 × (2 ** 10) = 131072
   
   # Exact in Python integers ✓

**Issue with Non-Integer Factors**:

.. code-block:: text

   nx = 128
   scale = 1.5
   
   # After 1 step
   nx_1 = int(128 × 1.5) = 192  # Truncation
   
   # After 2 steps (two separate scalings)
   nx_2 = int(192 × 1.5) = 288  # OK
   
   # Direct calculation
   nx_2_direct = int(128 × 1.5^2) = int(288.0) = 288 ✓

**Recommendation**: Use integer scale factors (2, 4, 8) for exact results.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PART 5: TESTING & VALIDATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

5.1 Unit Tests
==============

Test: Constant Work Per Process
--------------------------------

.. code-block:: python
   :linenos:

   def test_weak_scaling_constant_work():
       """Verify weak scaling maintains constant work per process."""
       
       # Setup
       engine = ScalingEngine(
           input_file="test_input.inp",
           output_dir="/tmp/test_output",
           nodes=[1, 2, 4, 8, 16],
           procs_per_node=8,
           scale_factor=2,
           dims=3,
           initial_procs=(2, 2, 2),
           initial_cells=(128, 128, 128),
           scaling_type="weak"
       )
       
       # Generate configurations
       configs = engine.generate_scaling_configs()
       
       # Calculate baseline cells per process
       baseline = configs[0]
       baseline_cells_per_proc = (
           baseline['nx'] * baseline['ny'] * baseline['nz']
       ) / (
           baseline['nproc_x'] * baseline['nproc_y'] * baseline['nproc_z']
       )
       
       # Verify all configurations have same work
       for config in configs:
           total_cells = config['nx'] * config['ny'] * config['nz']
           total_procs = (
               config['nproc_x'] * 
               config['nproc_y'] * 
               config['nproc_z']
           )
           cells_per_proc = total_cells / total_procs
           
           assert abs(cells_per_proc - baseline_cells_per_proc) < 0.01, \
               f"Node {config['nodes']}: " \
               f"cells_per_proc={cells_per_proc} != " \
               f"baseline={baseline_cells_per_proc}"

Test: Constant Grid Resolution
-------------------------------

.. code-block:: python
   :linenos:

   def test_weak_scaling_constant_resolution():
       """Verify grid resolution stays constant."""
       
       engine = ScalingEngine(
           input_file="test_input.inp",
           output_dir="/tmp/test_output",
           nodes=[1, 2, 4, 8],
           procs_per_node=8,
           scale_factor=2,
           dims=3,
           initial_procs=(2, 2, 2),
           initial_domain=(10.0, 10.0, 10.0),
           initial_cells=(128, 128, 128),
           scaling_type="weak"
       )
       
       configs = engine.generate_scaling_configs()
       
       # Calculate baseline resolution
       baseline = configs[0]
       dx_baseline = baseline['Lx'] / baseline['nx']
       dy_baseline = baseline['Ly'] / baseline['ny']
       dz_baseline = baseline['Lz'] / baseline['nz']
       
       # Verify all configurations
       for config in configs:
           dx = config['Lx'] / config['nx']
           dy = config['Ly'] / config['ny']
           dz = config['Lz'] / config['nz']
           
           assert abs(dx - dx_baseline) < 1e-10, \
               f"Node {config['nodes']}: dx changed"
           assert abs(dy - dy_baseline) < 1e-10, \
               f"Node {config['nodes']}: dy changed"
           assert abs(dz - dz_baseline) < 1e-10, \
               f"Node {config['nodes']}: dz changed"

Test: Dimension Rotation Pattern
---------------------------------

.. code-block:: text
   :linenos:

   def test_3d_dimension_rotation():
       """Verify correct X→Y→Z rotation in 3D."""
       
       engine = ScalingEngine(
           input_file="test_input.inp",
           output_dir="/tmp/test_output",
           nodes=[1, 2, 4, 8, 16, 32, 64],  # 7 steps
           procs_per_node=8,
           scale_factor=2,
           dims=3,
           initial_procs=(2, 2, 2),
           initial_cells=(128, 128, 128),
           scaling_type="weak"
       )
       
       configs = engine.generate_scaling_configs()
       
       # Expected pattern: X, Y, Z, X, Y, Z
       expected_topologies = [
           (2, 2, 2),  # Node 1: baseline
           (4, 2, 2),  # Node 2: X scaled
           (4, 4, 2),  # Node 4: Y scaled
           (4, 4, 4),  # Node 8: Z scaled
           (8, 4, 4),  # Node 16: X scaled again
           (8, 8, 4),  # Node 32: Y scaled again
           (8, 8, 8),  # Node 64: Z scaled again
       ]
       
       for i, config in enumerate(configs):
           actual = (
               config['nproc_x'],
               config['nproc_y'],
               config['nproc_z']
           )
           expected = expected_topologies[i]
           
           assert actual == expected, \
               f"Node {config['nodes']}: " \
               f"topology={actual} != expected={expected}"

5.2 Integration Tests
=====================

End-to-End Test
---------------

.. code-block:: python
   :linenos:

   def test_complete_weak_scaling_workflow():
       """Test complete weak scaling workflow."""
       
       import tempfile
       import os
       
       # Create temporary directories
       with tempfile.TemporaryDirectory() as tmpdir:
           input_dir = Path(tmpdir) / "input"
           output_dir = Path(tmpdir) / "output"
           input_dir.mkdir()
           output_dir.mkdir()
           
           # Create test input file
           input_file = input_dir / "test.inp"
           input_file.write_text("""
           # Test input file
           Lx = 10.0
           Ly = 10.0
           Lz = 10.0
           nxc = 128
           nyc = 128
           nzc = 128
           XLEN = 2
           YLEN = 2
           ZLEN = 2
           """)
           
           # Run scaling engine
           engine = ScalingEngine(
               input_file=str(input_file),
               output_dir=str(output_dir),
               nodes=[1, 2, 4],
               procs_per_node=8,
               scale_factor=2,
               dims=3,
               initial_procs=(2, 2, 2),
               initial_domain=(10.0, 10.0, 10.0),
               initial_cells=(128, 128, 128),
               scaling_type="weak"
           )
           
           configs = engine.generate_scaling_configs()
           
           # Verify configurations generated
           assert len(configs) == 3
           
           # Verify input files created
           for node_count in [1, 2, 4]:
               input_file = output_dir / f"input_n{node_count}.inp"
               assert input_file.exists(), \
                   f"Input file not created for node {node_count}"
               
               # Verify content
               content = input_file.read_text()
               assert "Lx" in content
               assert "XLEN" in content
           
           # Verify scaling correctness
           # ... (add assertions from previous tests)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PART 6: PRACTICAL USAGE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

6.1 Command-Line Usage
======================

Basic 3D Weak Scaling
----------------------

.. code-block:: bash

   python hpc_auto.py /path/to/code \\
       --scaling weak \\
       --nodes 64 \\
       --scaling-factor 2 \\
       --scaling-dimensions 3 \\
       --initial-procs 2,2,2 \\
       --initial-domain 10.0,10.0,10.0 \\
       --initial-cells 128,128,128

2D Weak Scaling (Fixed Z)
--------------------------

.. code-block:: bash

   python hpc_auto.py /path/to/code \\
       --scaling weak \\
       --scaling-dimensions 2 \\
       --initial-procs 2,2,2 \\
       --nodes 32

6.2 Python API Usage
====================

Detailed Example
----------------

.. code-block:: text
   :linenos:

   from pathlib import Path
   from engine.scaling import ScalingEngine
   from utils.report_generator import ReportGenerator
   
   # Configure scaling engine
   engine = ScalingEngine(
       input_file="simulation/input.inp",
       output_dir="output/weak_scaling",
       nodes=[1, 2, 4, 8, 16, 32, 64],  # Powers of 2
       procs_per_node=128,               # Processes per node
       scale_factor=2,                   # Double each step
       dims=3,                           # 3D scaling
       
       # Initial configuration (Node 1)
       initial_procs=(4, 4, 4),          # 64 total processes
       initial_domain=(50.0, 50.0, 50.0),
       initial_cells=(512, 512, 512),
       particles_per_cell=(100, 100, 100),
       
       scaling_type="weak"
   )
   
   # Generate scaling configurations
   print("Generating weak scaling configurations...")
   configs = engine.generate_scaling_configs()
   
   # Print configuration table
   print("\\n" + "="*80)
   print("Weak Scaling Configuration Table")
   print("="*80)
   print(f"{'Nodes':>6} {'Procs':>8} {'MPI Topology':>20} "
         f"{'Domain Size':>25} {'Grid Cells':>20}")
   print("-"*80)
   
   for config in configs:
       print(f"{config['nodes']:>6} "
             f"{config['total_procs']:>8} "
             f"{config['nproc_x']:>6}×{config['nproc_y']:<6}×{config['nproc_z']:<6} "
             f"{config['Lx']:>7.2f}×{config['Ly']:<7.2f}×{config['Lz']:<7.2f} "
             f"{config['nx']:>6}×{config['ny']:<6}×{config['nz']:<6}")
   
   print("="*80)
   
   # Verify work is constant
   print("\\nVerifying constant work per process...")
   baseline_work = (
       configs[0]['nx'] * configs[0]['ny'] * configs[0]['nz']
   ) / (
       configs[0]['nproc_x'] * 
       configs[0]['nproc_y'] * 
       configs[0]['nproc_z']
   )
   
   for config in configs:
       work = (config['nx'] * config['ny'] * config['nz']) / \\
              (config['nproc_x'] * config['nproc_y'] * config['nproc_z'])
       
       if abs(work - baseline_work) > 0.01:
           print(f"  ⚠ Node {config['nodes']}: "
                 f"work mismatch ({work} != {baseline_work})")
       else:
           print(f"  ✓ Node {config['nodes']}: constant work verified")

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
APPENDIX A: MATHEMATICAL PROOFS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Proof: Constant Work Per Process
=================================

**Theorem**: The dimension rotation algorithm maintains constant cells per process.

**Proof by Induction**:

*Base Case* (i = 0, Node 1):
  Let baseline cells per process = C₀
  
  C₀ = (N_x^0 × N_y^0 × N_z^0) / (P_x^0 × P_y^0 × P_z^0)

*Inductive Step*:
  Assume C_{i-1} = C₀
  
  At step i, one dimension d is scaled by factor s:
  - P_d^i = P_d^{i-1} × s
  - N_d^i = N_d^{i-1} × s
  - All other dimensions unchanged
  
  Therefore:
  
  C_i = (N_x^i × N_y^i × N_z^i) / (P_x^i × P_y^i × P_z^i)
      = (N_d^{i-1}×s × ∏_{j≠d} N_j^{i-1}) / (P_d^{i-1}×s × ∏_{j≠d} P_j^{i-1})
      = (N_d^{i-1} × ∏_{j≠d} N_j^{i-1}) / (P_d^{i-1} × ∏_{j≠d} P_j^{i-1})  [s cancels]
      = C_{i-1}
      = C₀  [by inductive hypothesis]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
REFERENCES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[1] Source Code: ``engine/scaling.py`` - Complete implementation
[2] Test Suite: ``tests/test_scaling.py`` - Verification tests
[3] User Guide: :doc:`../user_guide` - High-level usage
[4] API Reference: :doc:`../api/engine` - API documentation

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
END OF DOCUMENT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
