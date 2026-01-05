import math

Lx = 512.0
Ly = 128.0
n0 = 1.0
Tsim = 600.0
resx = 10.
resy = 10.
dx = 1/resx
dy = 1/resy
t0 = 0.95*dx/math.sqrt(2)
pusher = "norm"

def n0_(x,y):
        if (x>60):
                return 0.5*n0
        else:
                return 0.

def my_filter(x, y, px, py, pz):
    return (16.5<x)*(x<17.5)*(63<y)*(y<65)

Main(
    geometry = "2d3v",
    
    interpolation_order = 2 ,
    
    cell_length = [dx, dy],
    sim_length  = [Lx, Ly],
    
    number_of_patches = [8, 8],
    
    timestep = t0,
    sim_time = Tsim,
    
    bc_em_type_x = ['silver-muller','silver-muller'],
    bc_em_type_y = ['periodic', 'periodic'],

    print_every = 1,    

    random_seed = 0
)


ExtField(
    field = "Bz",
    profile = constant(1.0, xvacuum=60.0)
)

Species(
    species_type = "beamelectron",
    initPosition_type = "regular",
    initMomentum_type = "maxwell-juettner",
    n_part_per_cell = 16,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    nb_density = gaussian(n0, xvacuum=2., xlength=30., xfwhm=15., xcenter=17., xorder=2., yvacuum=32, ylength=64., yfwhm=32., ycenter=64., yorder=2.),
    mean_velocity = [0.9998, 0.0, 0.0],
    temperature = [1e-2],
    dynamics_type = pusher,
    bc_part_type_xmin  = "none",
    bc_part_type_xmax  = "none",
    bc_part_type_ymin = "none",
    bc_part_type_ymax = "none",
    isTest = False,
    track_every = 1,
    track_flush_every = 25,
    track_filter = my_filter,
)

Species(
    species_type = "beampositron",
    initPosition_type = "regular",
    initMomentum_type = "maxwell-juettner",
    n_part_per_cell = 16,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    nb_density = gaussian(n0, xvacuum=2., xlength=30., xfwhm=15., xcenter=17., xorder=2., yvacuum=32, ylength=64., yfwhm=32., ycenter=64., yorder=2.),
    mean_velocity = [0.9998, 0.0, 0.0],
    temperature = [1e-2],
    dynamics_type = pusher,
    bc_part_type_xmin  = "none",
    bc_part_type_xmax  = "none",
    bc_part_type_ymin = "none",
    bc_part_type_ymax = "none",
    isTest = False,
    track_every = 1,
    track_flush_every = 25,
    track_filter = my_filter,

)

Species(
    species_type = "plasmaelectron",
    initPosition_type = "regular",
    initMomentum_type = "maxwell-juettner",
    n_part_per_cell = 16,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    nb_density = n0_,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [1e-2],
    dynamics_type = pusher,
    bc_part_type_xmin  = "none",
    bc_part_type_xmax  = "none",
    bc_part_type_ymin = "none",
    bc_part_type_ymax = "none",
    isTest = False,
)

Species(
    species_type = "plasmaproton",
    initPosition_type = "regular",
    initMomentum_type = "maxwell-juettner",
    n_part_per_cell = 16,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    nb_density = n0_,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [1e-2],
    dynamics_type = pusher,
    bc_part_type_xmin  = "none",
    bc_part_type_xmax  = "none",
    bc_part_type_ymin = "none",
    bc_part_type_ymax = "none",
    isTest = False,
)


DiagFields(
    every = 250,
    fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz","Rho_beamelectron", "Rho_beampositron" ,"Rho_plasmaelectron", "Rho_plasmaproton"]
)

DiagScalar(
    every = 250 ,
    vars = ["Utot", "Ukin", "Ukin_beamelectron", "Ukin_beampositron", "Ukin_plasmaelectron", "Ukin_plasmaproton"],
    precision = 10,
)




