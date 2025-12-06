================================================================
Domain Decomposition: Complete Technical Reference
================================================================

.. contents:: Table of Contents
   :local:
   :depth: 3

Introduction for Everyone
===========================

What is Domain Decomposition?
------------------------------

**Domain decomposition** is the method of splitting a large problem into smaller pieces 
so multiple computers can work on it simultaneously.

**Simple Analogy: Painting a House**

.. code-block:: text

   You need to paint a large house.
   
   Option 1: One Painter
   - Paints the entire house alone
   - Takes 10 days
   
   Option 2: Domain Decomposition
   - Divide house into 4 sections (domain decomposition!)
   - Assign 1 painter to each section
   - Each works independently
   - Done in 2.5 days (4× faster!)
   
   The "domain" is the house.
   "Decomposition" is splitting it into sections.

Real Scientific Example
~~~~~~~~~~~~~~~~~~~~~~~~

**Weather Simulation Over North America**

.. code-block:: text

   Domain: North America (3000 km × 2000 km × 20 km atmosphere)
   
   Problem: Too large for one computer!
   - 10 billion grid points
   - Need 1TB of memory
   - Would take 1 year on 1 computer
   
   Solution: Domain Decomposition
   - Split into 1000 smaller regions
   - Each computer handles one region
   - Each needs only 1GB memory
   - Done in ~1 day using 1000 computers!

Why is This Important?
-----------------------

**Three Key Benefits**:

1. **Solve Bigger Problems**
   - One computer: Limited memory (64 GB typical)
   - 100 computers: 6,400 GB combined!
   - Can simulate much larger domains

2. **Get Results Faster**
   - Parallel processing
   - 100 computers can work simultaneously
   - Near-100× speedup (if done correctly!)

3. **Use Available Resources**
   - Most supercomputers have 1000+ nodes
   - Domain decomposition lets you use them all
   - Get your money's worth!

Conceptual Foundation
======================

Types of Decomposition
-----------------------

**1D Decomposition (Slicing)**

Like slicing bread:

.. code-block:: text

   Original domain: [============]
   
   Split into 4 slices:
   [===][===][===][===]
   
   Comp1 Comp2 Comp3 Comp4
   
   Good for: 1D problems (rivers, pipes, transmission lines)
   Limitation: Only 2 neighbors per process

**2D Decomposition (Checkerboard)**

Like a checkerboard:

.. code-block:: text

   Original domain: 
   [==========]
   [==========]
   [==========]
   [==========]
   
   Split into 4×4 grid:
   [==][==][==][==]
   [==][==][==][==]
   [==][==][==][==]
   [==][==][==][==]
   
   Good for: 2D problems (maps, surfaces, thin layers)
   Each process: Up to 4 neighbors

**3D Decomposition (Cubes)**

Like a Rubik's cube:

.. code-block:: text

   Original domain: Solid 3D cube
   
   Split into 2×2×2 = 8 smaller cubes:
   
   Front layer:         Back layer:
   [Cube1][Cube2]      [Cube5][Cube6]
   [Cube3][Cube4]      [Cube7][Cube8]
   
   Good for: 3D problems (atmosphere, ocean, solids)
   Each process: Up to 6 neighbors

Boundary Communication
-----------------------

**The Challenge: Processes Need to Talk**

.. code-block:: text

   When you split a domain, processes at boundaries need to share data:
   
   Example: Temperature simulation
   
   Process 1 domain:    Process 2 domain:
   [. . . . . 50°]      [55° . . . . .]
                ↕
         Boundary - needs communication!
   
   To calculate temperature at boundary:
   - Process 1 needs to know Process 2's value (55°)
   - Process 2 needs to know Process 1's value (50°)
   - They must exchange data!

**Ghost Cells / Halo Regions**

.. code-block:: text

   Solution: Each process keeps extra "ghost" cells
   
   Process 1's view:                Process 2's view:
   [. . . . . 50°][ghost]          [ghost][55° . . . . .]
                    ↑                 ↑
              Copy from P2        Copy from P1
   
   Before computation:
   1. Exchange boundary values
   2. Fill ghost cells
   3. Now each process can compute independently!

Load Balancing
---------------

**Equal Work for All Processes**

.. code-block:: text

   Bad decomposition:
   Process 1: 10,000 cells (takes 10 sec)
   Process 2: 100,000 cells (takes 100 sec)  ← Bottleneck!
   Process 3: 5,000 cells (takes 5 sec)
   
   Total time: 100 seconds (limited by slowest!)
   Processes 1 & 3 idle for 90+ seconds - wasted!
   
   Good decomposition:
   Process 1: 38,333 cells (takes 38.3 sec)
   Process 2: 38,333 cells (takes 38.3 sec)
   Process 3: 38,334 cells (takes 38.3 sec)
   
   Total time: 38.3 seconds
   All processes busy entire time!

Mathematical Foundations
=========================

Domain and Grid Definitions
-----------------------------

**Physical Domain**: :math:`\Omega \subset \mathbb{R}^3`

Represented by dimensions:

.. math::

   \Omega = [0, L_x] \times [0, L_y] \times [0, L_z]

Where:
- :math:`L_x, L_y, L_z` = domain lengths (e.g., kilometers)

**Computational Grid**: Discretization of :math:`\Omega`

.. math::

   n_x \times n_y \times n_z \text{ grid points}

**Grid Spacing**:

.. math::

   \Delta x = \frac{L_x}{n_x}, \quad \Delta y = \frac{L_y}{n_y}, \quad \Delta z = \frac{L_z}{n_z}

Decomposition Mathematics
---------------------------

**Given**:
- :math:`P` processes
- MPI topology :math:`(p_x, p_y, p_z)` where :math:`P = p_x \cdot p_y \cdot p_z`
- Grid :math:`(n_x, n_y, n_z)`

**Subdomain Size** for process at position :math:`(i_x, i_y, i_z)`:

.. math::

   n_x^{\text{local}} = \frac{n_x}{p_x}, \quad n_y^{\text{local}} = \frac{n_y}{p_y}, \quad n_z^{\text{local}} = \frac{n_z}{p_z}

**Subdomain Physical Size**:

.. math::

   L_x^{\text{local}} = \frac{L_x}{p_x}, \quad L_y^{\text{local}} = \frac{L_y}{p_y}, \quad L_z^{\text{local}} = \frac{L_z}{p_z}

**Grid Spacing** (same for all processes):

.. math::

   \Delta x^{\text{local}} = \Delta x = \frac{L_x}{n_x}

Implementation in HPC-ScaleTest
================================

Algorithm Overview
-------------------

.. code-block:: python

   class DomainDecomposer:
       """Decompose computational domain across MPI processes."""
       
       def __init__(self, domain_size, grid_dims, mpi_topology):
           """Initialize decomposer.
           
           Args:
               domain_size: (Lx, Ly, Lz) physical domain size
               grid_dims: (nx, ny, nz) grid dimensions
               mpi_topology: (px, py, pz) process grid
           """
           self.Lx, self.Ly, self.Lz = domain_size
           self.nx, self.ny, self.nz = grid_dims
           self.px, self.py, self.pz = mpi_topology
           
           # Verify divisibility
           assert self.nx % self.px == 0, "Grid not divisible in X"
           assert self.ny % self.py == 0, "Grid not divisible in Y"
           assert self.nz % self.pz == 0, "Grid not divisible in Z"
       
       def get_local_domain(self, process_id):
           """Get subdomain for a process.
           
           Returns:
               Dictionary with local domain information
           """
           # Get process coordinates
           ix, iy, iz = self._get_coords(process_id)
           
           # Calculate local grid size
           nx_local = self.nx // self.px
           ny_local = self.ny // self.py
           nz_local = self.nz // self.pz
           
           # Calculate local physical size
           Lx_local = self.Lx / self.px
           Ly_local = self.Ly / self.py
           Lz_local = self.Lz / self.pz
           
           # Calculate subdomain bounds
           x_start = ix * Lx_local
           x_end = (ix + 1) * Lx_local
           y_start = iy * Ly_local
           y_end = (iy + 1) * Ly_local
           z_start = iz * Lz_local
           z_end = (iz + 1) * Lz_local
           
           return {
               'process_id': process_id,
               'coordinates': (ix, iy, iz),
               'grid_size': (nx_local, ny_local, nz_local),
               'physical_size': (Lx_local, Ly_local, Lz_local),
               'bounds': {
                   'x': (x_start, x_end),
                   'y': (y_start, y_end),
                   'z': (z_start, z_end)
               }
           }

Worked Example
===============

Complete Domain Decomposition
-------------------------------

**Problem Setup**:

.. code-block:: text

   Simulating ocean currents in Pacific Ocean
   
   Physical domain:
   - Lx = 10,000 km (east-west)
   - Ly = 8,000 km (north-south)
   - Lz = 4 km (depth)
   
   Computational grid:
   - nx = 1000 points (10 km resolution in X)
   - ny = 800 points (10 km resolution in Y)
   - nz = 40 points (100 m resolution in Z)
   
   Total grid points: 1000 × 800 × 40 = 32,000,000
   
   Resources:
   - 64 processes available
   - MPI topology: (4, 4, 4) [calculated using MPI topology algorithm]

**Step 1: Calculate Local Grid Sizes**

.. code-block:: text

   Processes in each dimension:
   px = 4, py = 4, pz = 4
   
   Local grid size for each process:
   nx_local = 1000 / 4 = 250 points
   ny_local = 800 / 4 = 200 points
   nz_local = 40 / 4 = 10 points
   
   Each process handles: 250 × 200 × 10 = 500,000 points

**Step 2: Calculate Local Physical Domains**

.. code-block:: text

   Physical size for each process:
   Lx_local = 10,000 km / 4 = 2,500 km
   Ly_local = 8,000 km / 4 = 2,000 km
   Lz_local = 4 km / 4 = 1 km
   
   Each process simulates: 2,500 × 2,000 × 1 km region

**Step 3: Assign Subdomains to Processes**

.. code-block:: text

   Process 0 at (0,0,0):
   - X range: [0, 2500] km
   - Y range: [0, 2000] km
   - Z range: [0, 1] km
   - Grid: [0:250, 0:200, 0:10]
   
   Process 1 at (1,0,0):
   - X range: [2500, 5000] km
   - Y range: [0, 2000] km
   - Z range: [0, 1] km
   - Grid: [250:500, 0:200, 0:10]
   
   Process 17 at (1,0,1):
   - X range: [2500, 5000] km
   - Y range: [0, 2000] km
   - Z range: [1, 2] km
   - Grid: [250:500, 0:200, 10:20]
   
   [Continue for all 64 processes...]

**Verification**:

.. code-block:: text

   ✓ Each process has 500,000 points (equal distribution)
   ✓ All processes cover non-overlapping regions
   ✓ Together they cover entire domain:
     - X: 4 processes × 2,500 km = 10,000 km ✓
     - Y: 4 processes × 2,000 km = 8,000 km ✓
     - Z: 4 processes × 1 km = 4 km ✓

Summary
========

Key Takeaways
--------------

**Domain Decomposition**:
- Splits large problem into smaller subdomains
- Each process handles one subdomain
- Enables parallel computing

**Benefits**:
- Solve larger problems (distribute memory)
- Get faster results (parallel processing)
- Use available HPC resources efficiently

**HPC-ScaleTest Handles**:
- Automatic subdomain calculation
- Load balancing verification
- Boundary communication setup
- Configuration file generation

**Best Practices**:
- Choose process counts that evenly divide grid
- Balance subdomain sizes
- Minimize communication overhead
- Verify decomposition before running

This completes the domain decomposition technical reference!
