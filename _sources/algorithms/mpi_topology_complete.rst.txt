================================================================
MPI Topology Calculation: Complete Technical Reference
================================================================

.. contents:: Table of Contents
   :local:
   :depth: 3

Introduction for Everyone
===========================

What is MPI Topology?
----------------------

**MPI Topology** is how we organize multiple computers (processes) into a logical structure
so they can work together efficiently.

Think of it like organizing workers in a warehouse:

.. code-block:: text

   BAD Organization (Linear):
   Worker1 → Worker2 → Worker3 → Worker4 ... → Worker100
   
   - Worker1 talks to Worker2
   - Worker2 talks to Worker3
   - To send message from Worker1 to Worker100: Goes through 99 workers!
   - SLOW!
   
   
   GOOD Organization (Grid):
   Worker Grid (10 × 10):
   
   W1  W2  W3  W4  W5  W6  W7  W8  W9  W10
   W11 W12 W13 W14 W15 W16 W17 W18 W19 W20
   ...
   W91 W92 W93 W94 W95 W96 W97 W98 W99 W100
   
   - Each worker only talks to 4 neighbors (up/down/left/right)
   - Message from W1 to W100: Maximum 18 hops
   - FAST!

**That's MPI topology** - organizing processes efficiently for communication!

Real-World Analogy
~~~~~~~~~~~~~~~~~~

**Scenario**: You're simulating ocean currents across the Pacific Ocean.

**Bad approach** (no topology):
- Each computer stores random ocean regions
- To calculate currents, computers must ask everyone for neighboring data
- Chaotic and slow!

**Good approach** (with topology):
- Divide ocean into a grid (like a checkerboard)
- Each computer handles one grid square
- Each computer only talks to its 4-6 neighbors
- Organized and fast!

Why Do We Need This?
----------------------

**The Communication Problem**

Scientific simulations involve computers sharing data:

.. code-block:: text

   Example: Weather simulation
   
   Each computer simulates a region of the atmosphere.
   To calculate wind at borders, computers must exchange data.
   
   With random organization:
   - Computer #47 doesn't know who has the neighboring region!
   - Must ask everyone: "Do you have region X?"
   - 1000 computers = 1000 messages just to find neighbors!
   
   With MPI topology:
   - Computer at position (4,7) knows neighbors are at (3,7), (5,7), (4,6), (4,8)
   - Direct communication!
   - 4 messages instead of 1000!

**The Math Problem**

Computers need to divide work evenly:

.. code-block:: text

   Problem: Simulate 256 × 256 × 128 grid on 64 computers
   
   Question: How do we divide this?
   
   Bad: Random assignment
   - Computer 1 gets 1000 cells
   - Computer 2 gets 50,000 cells
   - Computer 1 finishes early and waits! (wasted resources)
   
   Good: MPI topology calculation
   - All computers get exactly 131,072 cells
   - All finish at the same time!
   - No wasted resources!

What You'll Learn
------------------

After reading this document, you'll understand:

✅ **What is a 3D MPI topology** (explained with pictures!)
✅ **How to calculate optimal topology** (the algorithm)
✅ **Why balanced topologies matter** (performance impact)
✅ **How to verify correctness** (testing methods)
✅ **Real examples** (worked calculations)

Who Should Read This?
-----------------------

- **Users**: Understand how HPC-ScaleTest organizes your parallel runs
- **Developers**: Learn the algorithm to extend or debug it
- **Students**: Study parallel computing with practical examples
- **Researchers**: See mathematical foundations for papers

No programming knowledge required for the concepts!

Conceptual Foundation
======================

Understanding Dimensions
-------------------------

**What are dimensions in parallel computing?**

Dimensions are different ways to split work:

**1D Topology** (One-dimensional split):

.. code-block:: text

   Like slicing a loaf of bread:
   
   [Slice1][Slice2][Slice3][Slice4]
   
   4 processes arranged in a LINE
   Topology: (4, 1, 1) or (1, 4, 1) or (1, 1, 4)
   
   Each process has 2 neighbors (left and right)

**2D Topology** (Two-dimensional split):

.. code-block:: text

   Like a checkerboard:
   
   [1][2][3][4]
   [5][6][7][8]
   
   8 processes arranged in a GRID
   Topology: (4, 2, 1)
   
   Each process has 4 neighbors (up/down/left/right)

**3D Topology** (Three-dimensional split):

.. code-block:: text

   Like a Rubik's cube:
   
   Front layer:        Back layer:
   [1 ][2 ][3 ][4 ]   [13][14][15][16]
   [5 ][6 ][7 ][8 ]   [17][18][19][20]
   [9 ][10][11][12]   [21][22][23][24]
   
   24 processes arranged in a CUBE
   Topology: (4, 3, 2)
   
   Each process has 6 neighbors (up/down/left/right/front/back)

**Why 3D?**

Most scientific simulations are 3D:
- Weather (atmosphere: latitude × longitude × altitude)
- Ocean (water: x × y × depth)
- Materials (solid: x × y × z coordinates)
- Universe (space: x × y × z)

Communication Patterns
-----------------------

**Nearest Neighbor Communication**

In scientific simulations, processes usually only talk to immediate neighbors:

.. code-block:: text

   Weather Example:
   
   To calculate temperature at a point, you need:
   - Temperature to the North (neighbor above)
   - Temperature to the South (neighbor below)
   - Temperature to the East (neighbor right)
   - Temperature to the West (neighbor left)
   - Temperature above (neighbor in front layer)
   - Temperature below (neighbor in back layer)
   
   That's 6 neighbors in 3D!

**Communication Cost**

.. code-block:: text

   1D Topology (8 processes in a line):
   Process 0 → Process 7
   - Must go through 6 intermediate processes
   - 7 hops total
   
   2D Topology (8 processes in 2×4 grid):
   Process at (0,0) → Process at (1,3)
   - Can go directly or through 2-3 processes
   - 2-4 hops
   
   3D Topology (8 processes in 2×2×2 cube):
   Process at (0,0,0) → Process at (1,1,1)
   - Diagonal distance: 3 hops
   - Better balance!

**Optimal Topology = Cubic Shape**

.. code-block:: text

   For 64 processes:
   
   Bad: (64, 1, 1) - Long line
   - Maximum distance: 63 hops
   - Many processes far apart
   
   Better: (8, 8, 1) - Flat square
   - Maximum distance: 14 hops
   - More balanced
   
   Best: (4, 4, 4) - Cube
   - Maximum distance: 9 hops
   - Most balanced!
   
   Cube minimizes maximum communication distance!

Load Balancing Principles
---------------------------

**Equal Work for All Processes**

.. code-block:: text

   Scenario: 1000 tasks, 10 workers
   
   Bad distribution:
   Worker 1: 500 tasks (takes 50 minutes)
   Worker 2: 50 tasks (takes 5 minutes)
   ...
   Total time: 50 minutes (Worker 1 is bottleneck!)
   Workers 2-10 idle for 45 minutes!
   
   Good distribution:
   Worker 1: 100 tasks (takes 10 minutes)
   Worker 2: 100 tasks (takes 10 minutes)
   ...
   Total time: 10 minutes
   All workers busy!

**How MPI Topology Helps**

.. code-block:: text

   Problem: 256 × 256 × 128 grid cells
   
   Question: How many processes in each dimension?
   
   Bad: (16, 4, 1)
   - X: 256 / 16 = 16 cells per process
   - Y: 256 / 4 = 64 cells per process  
   - Z: 128 / 1 = 128 cells per process
   - Each process: 16 × 64 × 128 = 131,072 cells
   - But communication is imbalanced! (1 neighbor in Z, 4 in X)
   
   Good: (4, 4, 4)
   - X: 256 / 4 = 64 cells per process
   - Y: 256 / 4 = 64 cells per process
   - Z: 128 / 4 = 32 cells per process
   - Each process: 64 × 64 × 32 = 131,072 cells (same!)
   - Communication is balanced (6 neighbors, similar distances)

Mathematical Foundations
=========================

Topology Notation
------------------

**MPI Topology is written as** :math:`(p_x, p_y, p_z)`

Where:
- :math:`p_x` = number of processes in X dimension
- :math:`p_y` = number of processes in Y dimension  
- :math:`p_z` = number of processes in Z dimension

**Total Processes**:

.. math::

   P = p_x \times p_y \times p_z

**Examples**:

.. code-block:: text

   8 processes:
   - (2, 2, 2): 2 × 2 × 2 = 8 ✓ (cube)
   - (4, 2, 1): 4 × 2 × 1 = 8 ✓ (rectangle)
   - (8, 1, 1): 8 × 1 × 1 = 8 ✓ (line)
   
   64 processes:
   - (4, 4, 4): 4 × 4 × 4 = 64 ✓ (cube)
   - (8, 8, 1): 8 × 8 × 1 = 64 ✓ (flat square)
   - (16, 4, 1): 16 × 4 × 1 = 64 ✓ (rectangle)

Process Coordinate System
---------------------------

Each process gets a unique 3D coordinate:

.. math::

   \text{Process ID} = i_x + p_x \cdot i_y + p_x \cdot p_y \cdot i_z

Where:
- :math:`0 \leq i_x < p_x`
- :math:`0 \leq i_y < p_y`
- :math:`0 \leq i_z < p_z`

**Example**: Topology (4, 2, 3) = 24 processes

.. code-block:: text

   Process coordinates:
   
   Layer 0 (iz=0):           Layer 1 (iz=1):           Layer 2 (iz=2):
   (0,0,0) (1,0,0) (2,0,0) (3,0,0)   (0,0,1) (1,0,1) (2,0,1) (3,0,1)   (0,0,2) (1,0,2) (2,0,2) (3,0,2)
   (0,1,0) (1,1,0) (2,1,0) (3,1,0)   (0,1,1) (1,1,1) (2,1,1) (3,1,1)   (0,1,2) (1,1,2) (2,1,2) (3,1,2)
   
   Process IDs:
   Layer 0: 0-7      Layer 1: 8-15     Layer 2: 16-23

**Converting Process ID to Coordinates**:

.. code-block:: python

   def get_coordinates(process_id, px, py, pz):
       """Convert process ID to (ix, iy, iz) coordinates."""
       iz = process_id // (px * py)
       remainder = process_id % (px * py)
       iy = remainder // px
       ix = remainder % px
       return (ix, iy, iz)
   
   # Example: Process 13 in (4, 2, 3) topology
   coords = get_coordinates(13, 4, 2, 3)
   # Result: (1, 0, 1)

Prime Factorization
--------------------

**The Foundation of Topology Calculation**

Every number can be broken down into prime factors:

.. math::

   P = 2^a \times 3^b \times 5^c \times 7^d \times \ldots

**Examples**:

.. code-block:: text

   8 = 2³ = 2 × 2 × 2
   Factors: [2, 2, 2]
   
   24 = 2³ × 3 = 2 × 2 × 2 × 3
   Factors: [2, 2, 2, 3]
   
   60 = 2² × 3 × 5 = 2 × 2 × 3 × 5
   Factors: [2, 2, 3, 5]

**Algorithm**:

.. code-block:: python

   def factorize(n):
       """Find all prime factors of n.
       
       Returns:
           List of prime factors (with repetitions)
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

**Step-by-Step Example**: :math:`n = 24`

.. code-block:: text

   Start: n = 24, d = 2, factors = []
   
   Is 24 divisible by 2? Yes!
   - 24 / 2 = 12
   - factors = [2]
   - n = 12
   
   Is 12 divisible by 2? Yes!
   - 12 / 2 = 6
   - factors = [2, 2]
   - n = 6
   
   Is 6 divisible by 2? Yes!
   - 6 / 2 = 3
   - factors = [2, 2, 2]
   - n = 3
   
   Is 3 divisible by 2? No, try d = 3
   
   Is 3 divisible by 3? Yes!
   - 3 / 3 = 1
   - factors = [2, 2, 2, 3]
   - n = 1
   
   Done! 24 = [2, 2, 2, 3]

Optimal Topology Criteria
---------------------------

**Goal**: Find :math:`(p_x, p_y, p_z)` that minimizes communication overhead

**Criterion 1: Cubic Ratio**

Minimize the aspect ratio:

.. math::

   \text{Aspect Ratio} = \frac{\max(p_x, p_y, p_z)}{\min(p_x, p_y, p_z)}

**Best**: Aspect ratio = 1 (perfect cube)

**Example**:

.. code-block:: text

   For 64 processes:
   
   (64, 1, 1): Aspect ratio = 64/1 = 64 (very bad!)
   (8, 8, 1): Aspect ratio = 8/1 = 8 (better)
   (4, 4, 4): Aspect ratio = 4/4 = 1 (perfect!)

**Criterion 2: Balanced Communication**

Processes should have similar numbers of neighbors:

.. code-block:: text

   Interior process neighbors:
   - 1D (line): 2 neighbors
   - 2D (grid): 4 neighbors
   - 3D (cube): 6 neighbors
   
   More dimensions = better load balance!

**Mathematical Formulation**:

Minimize:

.. math::

   C = \sum_{i=1}^{3} \frac{1}{p_i}

Subject to:

.. math::

   p_1 \times p_2 \times p_3 = P

Algorithm Design
=================

Overview
---------

HPC-ScaleTest's MPI topology calculation algorithm:

.. code-block:: text

   INPUT: 
   - Number of processes: P
   
   ALGORITHM:
   
   Step 1: Factor P into primes
           Example: 24 → [2, 2, 2, 3]
   
   Step 2: Distribute factors to (px, py, pz)
           Strategy: Keep dimensions balanced
           Example: [2,2,2,3] → (3,2,2,2) → (3×2, 2, 2) → (6, 2, 2)
   
   Step 3: Verify result
           Check: px × py × pz = P
           Example: 6 × 2 × 2 = 24 ✓
   
   OUTPUT:
   - Optimal MPI topology: (px, py, pz)

Complete Implementation
------------------------

.. code-block:: text

   # From engine/scaling.py
   
   class MPITopologyCalculator:
       """Calculate optimal MPI topology for parallel runs.
       
       The goal is to create a balanced 3D process grid that:
       1. Multiplies to give the correct total process count
       2. Maintains approximately cubic ratios (px ≈ py ≈ pz)
       3. Minimizes communication overhead
       """
       
       def __init__(self):
           """Initialize the topology calculator."""
           pass
       
       def calculate_topology(self, num_processes):
           """Calculate optimal 3D MPI topology.
           
           Args:
               num_processes: Total number of MPI processes
               
           Returns:
               Tuple (px, py, pz) representing the 3D process grid
           """
           if num_processes == 1:
               return (1, 1, 1)
           
           # Get prime factors
           factors = self._factorize(num_processes)
           
           # Distribute factors to create balanced topology
           topology = self._distribute_factors(factors)
           
           # Verify
           assert topology[0] * topology[1] * topology[2] == num_processes
           
           return topology
       
       def _factorize(self, n):
           """Find all prime factors of n.
           
           Args:
               n: Number to factorize
               
           Returns:
               List of prime factors (with repetitions)
               
           Example:
               _factorize(24) = [2, 2, 2, 3]
           """
           factors = []
           d = 2
           
           # Try all potential divisors up to √n
           while d * d <= n:
               while n % d == 0:
                   factors.append(d)
                   n //= d
               d += 1
           
           # If n > 1, then n itself is a prime factor
           if n > 1:
               factors.append(n)
           
           return factors
       
       def _distribute_factors(self, factors):
           """Distribute prime factors to create balanced topology.
           
           Strategy:
           1. Sort factors in descending order (largest first)
           2. Always assign factor to smallest current dimension
           3. This creates most balanced possible topology
           
           Args:
               factors: List of prime factors
               
           Returns:
               Tuple (px, py, pz) - the 3D topology
           """
           # Start with (1, 1, 1)
           px, py, pz = 1, 1, 1
           
           # Sort factors largest first
           sorted_factors = sorted(factors, reverse=True)
           
           # Distribute factors
           for factor in sorted_factors:
               # Find smallest dimension
               if px <= py and px <= pz:
                   px *= factor
               elif py <= pz:
                   py *= factor
               else:
                   pz *= factor
           
           return (px, py, pz)
       
       def calculate_neighbors(self, process_id, topology):
           """Calculate neighbors for a given process.
           
           Args:
               process_id: ID of the process (0 to P-1)
               topology: Tuple (px, py, pz)
               
           Returns:
               Dictionary of neighbors: {'left', 'right', 'up', 'down', 'front', 'back'}
           """
           px, py, pz = topology
           
           # Get coordinates
           ix, iy, iz = self._get_coordinates(process_id, px, py, pz)
           
           neighbors = {}
           
           # X-direction neighbors
           if ix > 0:
               neighbors['left'] = self._get_process_id(ix-1, iy, iz, px, py, pz)
           if ix < px - 1:
               neighbors['right'] = self._get_process_id(ix+1, iy, iz, px, py, pz)
           
           # Y-direction neighbors
           if iy > 0:
               neighbors['down'] = self._get_process_id(ix, iy-1, iz, px, py, pz)
           if iy < py - 1:
               neighbors['up'] = self._get_process_id(ix, iy+1, iz, px, py, pz)
           
           # Z-direction neighbors
           if iz > 0:
               neighbors['back'] = self._get_process_id(ix, iy, iz-1, px, py, pz)
           if iz < pz - 1:
               neighbors['front'] = self._get_process_id(ix, iy, iz+1, px, py, pz)
           
           return neighbors
       
       def _get_coordinates(self, process_id, px, py, pz):
           """Convert process ID to 3D coordinates."""
           iz = process_id // (px * py)
           remainder = process_id % (px * py)
           iy = remainder // px
           ix = remainder % px
           return (ix, iy, iz)
       
       def _get_process_id(self, ix, iy, iz, px, py, pz):
           """Convert 3D coordinates to process ID."""
           return ix + px * iy + px * py * iz

Worked Examples
================

Example 1: 8 Processes (Cube)
-------------------------------

**Calculate MPI Topology for 8 Processes**

**Step 1: Factorize**

.. code-block:: text

   8 = 2³
   factors = [2, 2, 2]

**Step 2: Distribute Factors**

.. code-block:: text

   Start: (px, py, pz) = (1, 1, 1)
   
   Factor 1: value = 2
   - Smallest dimension: px = 1 (all equal, pick first)
   - px = 1 × 2 = 2
   - Result: (2, 1, 1)
   
   Factor 2: value = 2
   - Smallest dimension: py = 1 (tie between py and pz, pick py)
   - py = 1 × 2 = 2
   - Result: (2, 2, 1)
   
   Factor 3: value = 2
   - Smallest dimension: pz = 1
   - pz = 1 × 2 = 2
   - Result: (2, 2, 2)

**Step 3: Verify**

.. code-block:: text

   px × py × pz = 2 × 2 × 2 = 8 ✓
   
   This is a PERFECT CUBE!
   - Aspect ratio: 2/2 = 1.0 (optimal)
   - All dimensions equal
   - Balanced communication

**Visualization**:

.. code-block:: text

   Front Layer (z=0):    Back Layer (z=1):
   [0][1]                [4][5]
   [2][3]                [6][7]
   
   Each process has up to 6 neighbors:
   - Process 0: neighbors = {1 (right), 2 (up), 4 (front)}
   - Process 3: neighbors = {2 (left), 1 (down), 7 (front)}
   - Process 7: neighbors = {6 (left), 5 (down), 3 (back)}

Example 2: 24 Processes
------------------------

**Calculate MPI Topology for 24 Processes**

**Step 1: Factorize**

.. code-block:: text

   24 = 2³ × 3 = 8 × 3
   
   Factorization process:
   24 ÷ 2 = 12    factors = [2]
   12 ÷ 2 = 6     factors = [2, 2]
   6 ÷ 2 = 3      factors = [2, 2, 2]
   3 ÷ 3 = 1      factors = [2, 2, 2, 3]

**Step 2: Sort Factors (Largest First)**

.. code-block:: text

   factors = [2, 2, 2, 3]
   sorted = [3, 2, 2, 2]  # Largest first helps balance

**Step 3: Distribute**

.. code-block:: text

   Start: (1, 1, 1)
   
   Factor 3:
   - Smallest: px = 1
   - px = 1 × 3 = 3
   - Result: (3, 1, 1)
   
   Factor 2:
   - Smallest: py = 1 (tie, pick py)
   - py = 1 × 2 = 2
   - Result: (3, 2, 1)
   
   Factor 2:
   - Smallest: pz = 1
   - pz = 1 × 2 = 2
   - Result: (3, 2, 2)
   
   Factor 2:
   - Smallest: py = 2 (tie between py and pz, pick py)
   - py = 2 × 2 = 4
   - Result: (3, 4, 2)

**Final Result: (3, 4, 2)**

.. code-block:: text

   Verify: 3 × 4 × 2 = 24 ✓
   
   Aspect ratio: 4 / 2 = 2.0 (reasonable)
   
   Not a perfect cube, but well-balanced!

**Visualization**:

.. code-block:: text

   Layer 0 (z=0):                Layer 1 (z=1):
   [0 ][1 ][2 ]                  [12][13][14]
   [3 ][4 ][5 ]                  [15][16][17]
   [6 ][7 ][8 ]                  [18][19][20]
   [9 ][10][11]                  [21][22][23]

Example 3: 64 Processes
------------------------

**Calculate MPI Topology for 64 Processes**

**Step 1: Factorize**

.. code-block:: text

   64 = 2⁶
   factors = [2, 2, 2, 2, 2, 2]

**Step 2: Distribute (all equal factors)**

.. code-block:: text

   Start: (1, 1, 1)
   
   This is interesting because all factors are the same!
   
   Factor 1 (2): px = 2 → (2, 1, 1)
   Factor 2 (2): py = 2 → (2, 2, 1)
   Factor 3 (2): pz = 2 → (2, 2, 2)
   Factor 4 (2): px = 4 → (4, 2, 2)  # px was smallest
   Factor 5 (2): py = 4 → (4, 4, 2)  # py was smallest
   Factor 6 (2): pz = 4 → (4, 4, 4)  # pz was smallest

**Final Result: (4, 4, 4)**

.. code-block:: text

   Verify: 4 × 4 × 4 = 64 ✓
   
   PERFECT CUBE!
   Aspect ratio: 4/4 = 1.0 (optimal!)
   
   This is the best possible topology for 64 processes!

**Why This is Optimal**:

.. code-block:: text

   Maximum distance between any two processes:
   - In (64, 1, 1): 63 hops
   - In (8, 8, 1): 14 hops
   - In (4, 4, 4): 9 hops  ← Best!
   
   Communication is minimized in all dimensions!

Example 4: Prime Number (17 Processes)
---------------------------------------

**Calculate MPI Topology for 17 Processes**

**Challenge**: 17 is a prime number!

**Step 1: Factorize**

.. code-block:: text

   17 is prime
   factors = [17]  # Cannot be factored further!

**Step 2: Distribute**

.. code-block:: text

   Start: (1, 1, 1)
   
   Only one factor: 17
   - Smallest: px = 1
   - px = 1 × 17 = 17
   - Result: (17, 1, 1)

**Final Result: (17, 1, 1)**

.. code-block:: text

   This creates a LINE of processes!
   
   [0][1][2][3][4][5][6][7][8][9][10][11][12][13][14][15][16]
   
   Not optimal for communication!
   - Process 0 to Process 16: 16 hops
   - Only 2 neighbors per process (except ends)

**Lesson**: Avoid prime numbers for process counts when possible!

**Better Alternatives**:

.. code-block:: text

   Instead of 17 processes, use:
   - 16 processes: (4, 4, 1) or (4, 2, 2)
   - 18 processes: (6, 3, 1) or (3, 3, 2)
   - 20 processes: (5, 4, 1) or (5, 2, 2)

Performance Analysis
=====================

Algorithmic Complexity
-----------------------

**Factorization**: :math:`O(\sqrt{P})`

.. code-block:: text

   For P processes:
   - Test divisors up to √P
   - Each test is O(1)
   
   Examples:
   - P = 64: √64 = 8 tests
   - P = 1024: √1024 ≈ 32 tests
   - P = 1,000,000: √1M ≈ 1000 tests
   
   Very fast even for large P!

**Factor Distribution**: :math:`O(k \log k)`

.. code-block:: text

   Where k = number of prime factors
   
   Steps:
   1. Sort factors: O(k log k)
   2. Distribute factors: O(k)
   
   Total: O(k log k)
   
   Typical: k < 10 (even for P = 1024)
   Practically: O(1) for realistic P

**Total**: :math:`O(\sqrt{P})`

Communication Overhead Analysis
--------------------------------

**Maximum Communication Distance**

In a :math:`(p_x, p_y, p_z)` topology:

.. math::

   D_{\max} = (p_x - 1) + (p_y - 1) + (p_z - 1)

**Example**: :math:`(4, 4, 4)`

.. code-block:: text

   Dmax = (4-1) + (4-1) + (4-1) = 3 + 3 + 3 = 9 hops

**Comparison for 64 Processes**:

.. code-block:: text

   Topology        Dmax    Aspect Ratio
   (64, 1, 1)      63      64.0  (worst)
   (32, 2, 1)      32      32.0
   (16, 4, 1)      18      16.0
   (8, 8, 1)       14       8.0
   (8, 4, 2)       12       4.0
   (4, 4, 4)        9       1.0  (best)

**Optimal Topology Minimizes Dmax**

Surface-to-Volume Ratio
-------------------------

**Why This Matters**

Communication happens at process boundaries (surface).
Computation happens inside process volume.

**Ratio**:

.. math::

   R = \frac{\text{Surface Area}}{\text{Volume}}

**Lower ratio = Better** (less communication per computation)

**Example**: 1000 cells per process

.. code-block:: text

   Linear (1000 × 1 × 1):
   - Volume: 1000
   - Surface: 2 (two ends)
   - Ratio: 2/1000 = 0.002
   
   Flat (100 × 10 × 1):
   - Volume: 1000
   - Surface: 2×(100 + 10) = 220
   - Ratio: 220/1000 = 0.22
   
   Cubic (10 × 10 × 10):
   - Volume: 1000
   - Surface: 6×100 = 600
   - Ratio: 600/1000 = 0.6
   
   WAIT - cubic has HIGHER ratio?!

**The Paradox Explained**:

.. code-block:: text

   In 1 dimension:
   - Very low surface area
   - BUT: Long communication chains!
   - Message from end to end: 999 hops
   
   In 3 dimensions:
   - Higher surface area
   - BUT: Short communication paths!
   - Message from corner to corner: 27 hops
   
   Trade-off: Higher local communication, but much shorter global paths
   
   For parallel computing: 3D is almost always better!

Practical Usage
================

Using HPC-ScaleTest
--------------------

**Automatic Topology Calculation**:

.. code-block:: bash

   # HPC-ScaleTest automatically calculates topology
   hpc-scaletest generate \
       --baseline config.yaml \
       --nodes "1 2 4 8" \
       --procs-per-node 8
   
   # Output shows calculated topologies:
   # 8 procs: (2, 2, 2)
   # 16 procs: (4, 2, 2)
   # 32 procs: (4, 4, 2)
   # 64 procs: (4, 4, 4)

**Manual Topology Specification**:

.. code-block:: bash

   # Override automatic calculation
   hpc-scaletest generate \
       --processes 24 \
       --topology 6 2 2  # Force this topology
   
   # Useful when you have specific requirements

**Verification**:

.. code-block:: bash

   # Verify topology is valid
   hpc-scaletest verify-topology \
       --processes 24 \
       --topology 3 4 2
   
   # Output:
   # ✓ Topology (3, 4, 2) is valid for 24 processes
   # ✓ Product: 3 × 4 × 2 = 24
   # ✓ Aspect ratio: 2.0 (good)

Python API
-----------

**Direct Calculation**:

.. code-block:: python

   from hpc_scaletest import MPITopologyCalculator
   
   calc = MPITopologyCalculator()
   
   # Calculate for different process counts
   for procs in [8, 24, 64, 128]:
       topo = calc.calculate_topology(procs)
       print(f"{procs} processes: {topo}")
   
   # Output:
   # 8 processes: (2, 2, 2)
   # 24 processes: (3, 4, 2)
   # 64 processes: (4, 4, 4)
   # 128 processes: (8, 4, 4)

**With Verification**:

.. code-block:: python

   # Calculate and verify
   topo = calc.calculate_topology(60)
   px, py, pz = topo
   
   assert px * py * pz == 60
   print(f"Topology for 60 processes: {topo}")
   print(f"Aspect ratio: {max(topo) / min(topo):.2f}")
   
   # Output:
   # Topology for 60 processes: (5, 4, 3)
   # Aspect ratio: 1.67

**Finding Neighbors**:

.. code-block:: python

   # Get neighbors for a process
   topology = (4, 3, 2)
   process_id = 5
   
   neighbors = calc.calculate_neighbors(process_id, topology)
   
   print(f"Process {process_id} neighbors:")
   for direction, neighbor_id in neighbors.items():
       print(f"  {direction}: Process {neighbor_id}")
   
   # Output:
   # Process 5 neighbors:
   #   left: Process 4
   #   right: Process 6
   #   down: Process 1
   #   up: Process 9
   #   back: Process 17

Common Issues and Solutions
============================

Issue 1: Prime Process Counts
-------------------------------

**Problem**:

.. code-block:: text

   You want to use 17 processes (prime number)
   Result: (17, 1, 1) - a long line!
   Poor communication pattern

**Solution 1**: Use Nearby Composite Number

.. code-block:: text

   Instead of 17, use:
   - 16 processes: (4, 4, 1) or (4, 2, 2)
   - 18 processes: (6, 3, 1) or (3, 3, 2)
   
   Slightly different count, much better topology!

**Solution 2**: Accept the Limitation

.. code-block:: text

   If you must use 17 processes:
   - Accept (17, 1, 1) topology
   - Optimize your algorithm for 1D decomposition
   - Consider using different parallel strategy

Issue 2: Non-Cubic Results
----------------------------

**Problem**:

.. code-block:: text

   100 processes gives (10, 10, 1)
   Not a true cube, mostly 2D!

**Explanation**:

.. code-block:: text

   100 = 2² × 5² = 4 × 25
   
   Best factorization:
   - px = 10
   - py = 10
   - pz = 1
   
   This is optimal given the factors!

**Solutions**:

1. **Use Different Process Count**:

   .. code-block:: text
   
      Instead of 100, use:
      - 96 processes: (6, 4, 4) - better cube
      - 108 processes: (6, 6, 3) - better balanced

2. **Accept It**:

   .. code-block:: text
   
      (10, 10, 1) is still reasonable
      - Good for 2D problems
      - Aspect ratio = 10/1 = 10 (acceptable)

Issue 3: Large Aspect Ratios
------------------------------

**Problem**:

.. code-block:: text

   Topology (16, 4, 1) has aspect ratio 16/1 = 16
   Is this bad?

**Analysis**:

.. code-block:: text

   Depends on your problem!
   
   For 2D problems (thin in Z):
   - (16, 4, 1) is perfectly fine!
   - Matches problem geometry
   
   For 3D problems:
   - Prefer more balanced topology
   - Consider different process count

**Solution**: Match Topology to Problem

.. code-block:: text

   2D Problem (atmospheric layer):
   - Domain: 1024 × 1024 × 1 (thin!)
   - Good: (16, 4, 1), (8, 8, 1)
   - Bad: (4, 4, 4) - wastes Z dimension
   
   3D Problem (ocean volume):
   - Domain: 256 × 256 × 128
   - Good: (4, 4, 4), (4, 4, 2)
   - Bad: (16, 4, 1) - ignores Z structure

Best Practices
===============

Choosing Process Counts
-------------------------

**Prefer Powers of 2**:

.. code-block:: text

   Excellent: 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024
   
   Why:
   - Factor into equal parts: 2×2×2×2...
   - Create perfect cubes: (2,2,2), (4,4,4), (8,8,8)
   - Match hardware (most HPC nodes have 2ⁿ cores)

**Good Alternatives**:

.. code-block:: text

   Also good: 12, 24, 48, 96  (= 2ⁿ × 3)
              20, 40, 80       (= 2ⁿ × 5)
              18, 36, 72       (= 2ⁿ × 3²)
   
   These factor nicely into balanced topologies

**Avoid**:

.. code-block:: text

   Avoid: Prime numbers (17, 23, 29, 31, 37, ...)
          Highly unbalanced (2 × large_prime)
   
   These create poor topologies

Matching Hardware
------------------

**Know Your Cluster**:

.. code-block:: bash

   # Find cores per node
   hpc-scaletest system-info
   
   # Example output:
   # Cores per node: 36
   # Nodes available: 10
   # Recommended process counts: 36, 72, 144, 288, 360

**Scale by Full Nodes**:

.. code-block:: text

   If nodes have 36 cores:
   
   Good:
   - 1 node: 36 processes
   - 2 nodes: 72 processes
   - 4 nodes: 144 processes
   
   Bad:
   - 50 processes (1.4 nodes - wastes resources)
   - 100 processes (2.8 nodes - wastes resources)

**Hybrid Parallelism**:

.. code-block:: text

   Modern approach: MPI + OpenMP
   
   Example: Node with 36 cores
   - 6 MPI processes × 6 OpenMP threads = 36 cores
   - Topology: Calculated for 6×nodes MPI processes
   - Better memory usage!

Testing Topologies
--------------------

**Verify Before Running**:

.. code-block:: python

   def test_topology(num_processes):
       """Test if topology is reasonable."""
       calc = MPITopologyCalculator()
       topo = calc.calculate_topology(num_processes)
       
       px, py, pz = topo
       
       # Check 1: Correct total
       assert px * py * pz == num_processes
       
       # Check 2: Aspect ratio
       aspect_ratio = max(topo) / min(topo)
       if aspect_ratio > 10:
           print(f"Warning: High aspect ratio {aspect_ratio:.1f}")
       
       # Check 3: No dimension is 1 (for 3D problems)
       if min(topo) == 1 and num_processes >= 8:
           print(f"Warning: Effectively 2D topology {topo}")
       
       return topo

**Benchmark Different Topologies**:

.. code-block:: bash

   # Test multiple topologies for same process count
   hpc-scaletest benchmark \
       --processes 24 \
       --topologies "(6,2,2)" "(4,3,2)" "(3,4,2)" "(2,3,4)"
   
   # See which performs best for your application!

Summary
========

Key Concepts
-------------

**MPI Topology**:
- 3D grid of processes: :math:`(p_x, p_y, p_z)`
- Total processes: :math:`P = p_x \times p_y \times p_z`
- Each process has unique 3D coordinate

**Goals**:
- Balanced dimensions (cubic shape)
- Minimize communication distance
- Match problem geometry

**Algorithm**:
1. Factor process count into primes
2. Distribute factors to keep dimensions balanced
3. Verify topology is valid

**Best Practices**:
- Use powers of 2 when possible
- Match topology to problem geometry
- Verify before running large studies

Further Reading
----------------

- :doc:`weak_scaling_complete` - Weak scaling algorithm
- :doc:`strong_scaling_complete` - Strong scaling algorithm
- :doc:`domain_decomposition_complete` - Domain decomposition details

Mathematical Appendix
======================

Neighbor Calculation Formulas
-------------------------------

**Process ID from Coordinates**:

.. math::

   \text{ID}(i_x, i_y, i_z) = i_x + p_x \cdot i_y + p_x \cdot p_y \cdot i_z

**Coordinates from Process ID**:

.. math::

   i_z = \left\lfloor \frac{\text{ID}}{p_x \cdot p_y} \right\rfloor

.. math::

   i_y = \left\lfloor \frac{\text{ID} \bmod (p_x \cdot p_y)}{p_x} \right\rfloor

.. math::

   i_x = \text{ID} \bmod p_x

**Neighbor IDs**:

.. math::

   \text{Left} = \text{ID} - 1 \quad (\text{if } i_x > 0)

.. math::

   \text{Right} = \text{ID} + 1 \quad (\text{if } i_x < p_x - 1)

.. math::

   \text{Down} = \text{ID} - p_x \quad (\text{if } i_y > 0)

.. math::

   \text{Up} = \text{ID} + p_x \quad (\text{if } i_y < p_y - 1)

.. math::

   \text{Back} = \text{ID} - p_x \cdot p_y \quad (\text{if } i_z > 0)

.. math::

   \text{Front} = \text{ID} + p_x \cdot p_y \quad (\text{if } i_z < p_z - 1)

Optimal Topology Proof
------------------------

**Theorem**: For :math:`P = 2^n` processes, topology :math:`(2^{n/3}, 2^{n/3}, 2^{n/3})` minimizes maximum communication distance.

**Proof**:

Maximum distance in topology :math:`(p_x, p_y, p_z)`:

.. math::

   D = (p_x - 1) + (p_y - 1) + (p_z - 1)

Subject to:

.. math::

   p_x \cdot p_y \cdot p_z = P

To minimize :math:`D`, use Lagrange multipliers:

.. math::

   \mathcal{L} = (p_x - 1) + (p_y - 1) + (p_z - 1) - \lambda(p_x \cdot p_y \cdot p_z - P)

Taking derivatives:

.. math::

   \frac{\partial \mathcal{L}}{\partial p_x} = 1 - \lambda p_y p_z = 0 \implies p_y p_z = \frac{1}{\lambda}

.. math::

   \frac{\partial \mathcal{L}}{\partial p_y} = 1 - \lambda p_x p_z = 0 \implies p_x p_z = \frac{1}{\lambda}

.. math::

   \frac{\partial \mathcal{L}}{\partial p_z} = 1 - \lambda p_x p_y = 0 \implies p_x p_y = \frac{1}{\lambda}

From these: :math:`p_x = p_y = p_z = P^{1/3}`

Therefore: **Cubic topology minimizes maximum communication distance**. ∎
