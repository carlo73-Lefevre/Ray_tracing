import os
import numpy as np
import sys
import pandas as pd
import scipy.io
import math
from datetime import datetime, timedelta
import re
from astropy.time import Time

import torch
import intel_extension_for_pytorch as ipex

device = torch.device('xpu' if torch.xpu.is_available() else 'cpu')
model = YourModel().to(device)

# Costanti
LIGHT_SPEED = 299792458  # Velocità della luce in m/s
SOLAR_CONSTANT = 1360.8  # Costante solare in W/m^2

# %% Funzioni Ratracing
def material_number(filename, path):
    return {"rho": 0.9}  # Valore di esempio

class SAT:
    def __init__(self, box_name, path):
        data = read_mesh_mat_txt(box_name, path)
        
        # Verifica e assegnazione dei dati caricati
        self.vertices = data.get("vertices", np.array([])) 
        self.faces = data.get("faces", np.array([]))
        self.material = data.get("material", "default_material")  
        self.normals_mesh = np.zeros_like(self.faces)  
        self.centers_mesh = np.zeros((self.faces.shape[0], 3)) if self.faces.size > 0 else np.array([])

    def rotate_panels(self, rotation_matrix):
        if self.vertices.size > 0:
            self.vertices = np.dot(self.vertices, rotation_matrix)
            self._calculate_normals()  

    def _calculate_normals(self):
        if self.faces.size > 0:
            for i, face in enumerate(self.faces):
                v1, v2, v3 = self.vertices[face]
                normal = np.cross(v2 - v1, v3 - v1)
                self.normals_mesh[i] = normal / np.linalg.norm(normal)  
            print("Calcolati i normali per ogni faccia.")

    def calculate_sun_angles(self, sun_vector):
        azimuth = np.arctan2(sun_vector[1], sun_vector[0])
        elevation = np.arcsin(sun_vector[2] / np.linalg.norm(sun_vector))
        return np.degrees(azimuth), np.degrees(elevation)

    def find_coplanar_faces(self):
        print("Identificazione delle facce complanari (metodo da definire).")


class c_Rays:
    def __init__(self):
        self.face_source = None  
        self.p_source = None     
        self.dir = None          
        self.face_from = None    
        self.p_from = None       
        self.dir_from = None     

def read_mesh_mat_txt(filename, path):
    data_path = os.path.join(path, "Dati Caomsol")
    os.chdir(data_path)

    S = {
        "vertices": np.array([]),
        "faces": np.array([]),
        "materials_idx": np.array([]),
        "optical_par": None,
        "Param": None,
        "material": "default_material"  
    }

    # Carica i dati dei materiali
    S["optical_par"] = material_number(filename, path)

    # Importa dati della mesh
    filename_mesh = f"{filename}.txt"
    if os.path.exists(filename_mesh):
        with open(filename_mesh, 'r') as file:
            lines = file.readlines()
            
        V = int(lines[4][20:].strip())
        F = int(lines[5][20:].strip())

        # Caricamento dati vertici e facce
        S["vertices"] = np.genfromtxt(lines[10:10+V], delimiter=' ', usecols=(0, 1, 2))
        S["faces"]    = np.genfromtxt(lines[10+V:10+V+F], delimiter=' ', usecols=(0, 1, 2)).astype(int)
        S["materials_idx"] = np.genfromtxt(lines[10+V+F:], delimiter=' ', usecols=0)

    filename_dat = f"{filename}.dat"
    if os.path.exists(filename_dat):
        opts = {"delim_whitespace": True, "header": None, "names": ["alfaA", "VarName2", "absorptionCoefficient"]}
        S["Param"] = pd.read_table(filename_dat, **opts)
        
    print(f"Dati caricati correttamente da {filename_mesh} e {filename_dat}")
    return S

def exclude_faces(Sat0, MAT, ray_directions, cont=1):
    IDX_Facce_ill = MAT["Dis_sort_ill"]

    for face_idx in IDX_Facce_ill:
        ray_dir = ray_directions[cont - 1]  
        face_vertices = Sat0.vertices[Sat0.faces[face_idx]]
        if not is_point_in_shadow(face_vertices, ray_dir):
            continue  
        else:
            IDX_Facce_ill = np.delete(IDX_Facce_ill, np.where(IDX_Facce_ill == face_idx))
    
    return IDX_Facce_ill

def is_point_in_shadow(face_vertices, ray_direction):
    return False

def calculate_srp_reflection(Sat, ray_directions, sun_dir, m_sat, FI_mod=SOLAR_CONSTANT):
    total_area = 0
    srp_forces = []

    for i, face in enumerate(Sat.faces):
        face_vertices = Sat.vertices[face]
        material_idx = Sat.materials_idx[i]
        rho = Sat.optical_par[material_idx].get("rho", 0.9)

        a, b, c = face_vertices
        AB = b - a
        AC = c - a
        area_face = 0.5 * np.linalg.norm(np.cross(AB, AC))
        total_area += area_face

        srp_force_face = area_face * rho * FI_mod / LIGHT_SPEED * sun_dir
        srp_forces.append(srp_force_face)

    return np.array(srp_forces)

# %% funzioni Visco

def sun(jdate):
    """
    Calcola l'ascensione retta, la declinazione e il vettore posizione del Sole in coordinate ECI.
    
    Args:
        jdate (float or np.array): Giorno giuliano.
        
    Returns:
        tuple: Ascensione retta (radians), declinazione (radians) e posizione ECI del Sole (km).
    """
    atr = np.pi / 648000
    rsun = np.zeros((3, len(np.atleast_1d(jdate))))

    # Calcolo delle variabili temporali
    djd = jdate - 2451545  # Giorni dal giorno giuliano di riferimento
    t = (djd / 36525) + 1  # Periodo in secoli giuliani

    # Argomenti fondamentali (in radianti)
    gs = r2r(0.993126 + 0.0027377785 * djd)
    lm = r2r(0.606434 + 0.03660110129 * djd)
    ls = r2r(0.779072 + 0.00273790931 * djd)
    g2 = r2r(0.140023 + 0.00445036173 * djd)
    g4 = r2r(0.053856 + 0.00145561327 * djd)
    g5 = r2r(0.056531 + 0.00023080893 * djd)
    rm = r2r(0.347343 - 0.00014709391 * djd)

    # Longitudine eclittica geocentrica del Sole
    plon = (6910 * np.sin(gs) + 72 * np.sin(2 * gs) - 17 * t * np.sin(gs) 
            - 7 * np.cos(gs - g5) + 6 * np.sin(lm - ls) 
            + 5 * np.sin(4 * gs - 8 * g4 + 3 * g5) 
            - 5 * np.cos(2 * (gs - g2)) - 4 * (np.sin(gs - g2) - np.cos(4 * gs - 8 * g4 + 3 * g5)) 
            + 3 * (np.sin(2 * (gs - g2)) - np.sin(g5) - np.sin(2 * (gs - g5))) )
    plon = ls + atr * (plon - 17 * np.sin(rm))

    # Distanza geocentrica del Sole in km
    rsm = 149597870.691 * (1.00014 - 0.01675 * np.cos(gs) - 0.00014 * np.cos(2 * gs))

    # Obliquità dell'eclittica
    obliq = atr * (84428 - 47 * t + 9 * np.cos(rm))

    # Ascensione retta e declinazione
    a = np.sin(plon) * np.cos(obliq)
    b = np.cos(plon)
    rasc = np.arctan2(a, b)
    decl = np.arcsin(np.sin(obliq) * np.sin(plon))

    # Vettore posizione ECI del Sole
    rsun[0, :] = rsm * np.cos(rasc) * np.cos(decl)
    rsun[1, :] = rsm * np.sin(rasc) * np.cos(decl)
    rsun[2, :] = rsm * np.sin(decl)

    return rasc, decl, rsun


def r2r(angle):
    """
    Converte un angolo in radianti nel range [0, 2*pi].

    Args:
        angle (float or np.array): Angolo in radianti.

    Returns:
        float or np.array: Angolo convertito nel range [0, 2*pi].
    """
    return angle % (2 * np.pi)


# Aggiorniamo la funzione principale per integrare questa nuova funzione "sun"
def Generazione_coordinate_Satellite_sole(satellite, dt, day_long):
    # Load constants and parameters
    par_directory = 'C:/Users/carlo/OneDrive - INAF - Istituto Nazionale di Astrofisica/G4S/Box_wing/Massimo_risultati/'
    
    const = read_par(par_directory + 'constant_SI.par')
    astrody = read_par(par_directory + 'astrody_SI.par')
    par_bw = read_box_wing_par(par_directory + 'box_wing_meta.par')


    # Satellite-specific parameters
    if satellite == 'E18':
        a = [27978099.6571506]  # Semi-major axis in meters
        e = [0.160352416634393]  # Eccentricity
        I = [50.369485965388 / 180 * np.pi]  # Inclination in radians
        OMG = [53.5051781754073 / 180 * np.pi]  # Longitude of ascending node in radians
        om = [50.1835975672884 / 180 * np.pi]  # Argument of perigee in radians
        MA = [316.069 / 180 * np.pi]  # Mean anomaly in radians
        dOMG = -0.04000414086 / 180 * np.pi  # RAAN rate in radians/day
        dom = 0.04910776939 / 180 * np.pi  # Argument of perigee rate in radians/day
        dMA = 667.909221051 / 180 * np.pi  # Mean anomaly rate in radians/day
        t_launch = mjuliandate([2015, 12, 17, 11, 51, 0]) # Launch date in MJD
        m_sat = 660.977  # Satellite mass in kg
    elif satellite == 'E08':
        a = [29600.351e3]  # Semi-major axis in meters
        e = [2.274e-4]  # Eccentricity
        I = [54.907 / 180 * np.pi]  # Inclination in radians
        OMG = [197.2821385041630 / 180 * np.pi]  # Longitude of ascending node in radians
        om = [360.0 / 180 * np.pi]  # Argument of perigee in radians
        MA = [130.77 / 180 * np.pi]  # Mean anomaly in radians
        dOMG = -0.02773706568296141 / 180 * np.pi  # RAAN rate in radians/day
        dom = 0.192 / 180 * np.pi  # Argument of perigee rate in radians/day
        dMA = 613.7 / 180 * np.pi  # Mean anomaly rate in radians/day
        t_launch = Time('2015-12-17 11:51:00', scale='utc').mjd  # Launch date in MJD
        m_sat = 719  # Satellite mass in kg

    # Initial Conditions Table
    Row_Names = [
        'semimajoraxes', 'eccentricity', 'Inclination',
        'longitude ascending node', 'pericenter argument',
        'mean anomaly', 'longitude ascending node variation',
        'pericenter argument rate', 'mean anomaly rate',
        'launch date', 'm_sat'
    ]
    init_Cond = pd.DataFrame({
        'VariableNames': Row_Names,
        'Values': [
            a[0], e[0], I[0], OMG[0], om[0], MA[0],
            dOMG, dom, dMA, t_launch, m_sat
        ]
    })

    # Time parameters
    t_ref = mjuliandate([2016,12,21, 00,00,00])  # Reference time in MJD
    t_start = t_launch  # Start time
    t_end = t_start + day_long  # End time
    t = np.arange(t_start, t_end + dt, dt)  # Time array
    t_par = t - t_ref  # Time difference from reference

    # Keplerian parameters over time
    OMG_t = OMG[0] + dOMG * t_par
    om_t = om[0] + dom * t_par
    MA_t = MA[0] + dMA * t_par
    n_orb = np.array([
        np.sin(I[0]) * np.sin(OMG_t),
        -np.sin(I[0]) * np.cos(OMG_t),
        np.cos(I[0]) * np.ones(len(OMG_t))
    ])

    # Satellite position in ECI frame
    r_sat = satellite_orb_ecc(a[0], om_t, MA_t, OMG_t, I[0], e[0])

    # Sun position in ECI frame
    t_jd = t + 2400000.5  # Convert MJD to JD
    times = Time(t_jd, format='jd', scale='utc')
    sun_coords = get_sun(times)
    sun_coords.representation_type = 'cartesian'
    sun_positions = np.array([
        sun_coords.x.to(u.m).value,
        sun_coords.y.to(u.m).value,
        sun_coords.z.to(u.m).value
    ])

    # Vector from satellite to sun
    r_sun_sat = sun_positions - r_sat

    # Distance calculations
    dr = np.linalg.norm(r_sun_sat, axis=0)  # Earth-Sun distance
    r_sun_norm = sun_positions / np.linalg.norm(sun_positions, axis=0)
    r_sat_m = np.linalg.norm(r_sat, axis=0)
    r_sat_norm = r_sat / r_sat_m
    r = r_sun_sat / dr

    # Shadow function
    nu = umbra(sun_positions, r_sat, astrody['RS'], astrody['RT'])

    # Satellite attitude in ECI frame
    Z_sat = -r_sat_norm  # Towards Earth center
    Y_sat = np.cross(r.T, Z_sat.T).T  # Perpendicular to sun position
    X_sat = np.cross(Y_sat.T, Z_sat.T).T
    X_sat = normalize_vectors(X_sat)
    Y_sat = normalize_vectors(Y_sat)
    Z_sat = normalize_vectors(Z_sat)

    X0 = np.cross(Z_sat.T, n_orb.T).T
    X0 = normalize_vectors(X0)

    # Yaw Steering Law
    Psi_nc = np.arctan2(
        np.einsum('ij,ij->j', r_sun_norm, n_orb),
        np.einsum('ij,ij->j', r_sun_norm, np.cross(r_sat_norm.T, n_orb.T).T)
    )

    # Auxiliary angles beta and epsilon
    beta = np.pi / 2 - angle_between_vectors(r_sun_norm, n_orb)
    x = np.cross(n_orb.T, r_sun_norm.T).T
    y = np.cross(n_orb.T, x.T).T
    eps = np.arccos(np.einsum('ij,ij->j', r_sat_norm, y))
    eps[eps > np.pi / 2] = np.pi - eps[eps > np.pi / 2]

    # Psi modification according to metadata
    Psi_meta, delta_Psi_meta = modify_psi(Psi_nc, beta, eps, dt)

    # Adjust satellite attitude due to modified Psi
    delta_Psi_meta = delta_Psi_meta.reshape(1, -1)
    X_sat_p = X_sat * np.cos(delta_Psi_meta) - Y_sat * np.sin(delta_Psi_meta)
    Y_sat_p = Y_sat * np.cos(delta_Psi_meta) + X_sat * np.sin(delta_Psi_meta)
    X_sat = X_sat_p
    Y_sat = Y_sat_p

    n_pan = np.cross(np.cross(Y_sat.T, r.T).T, Y_sat)
    n_pan = normalize_vectors(n_pan)

    # Solar vector in Body frame
    sun_body_x = np.einsum('ij,ij->j', r, X_sat)
    sun_body_y = np.einsum('ij,ij->j', r, Y_sat)
    sun_body_z = np.einsum('ij,ij->j', r, Z_sat)
    sun_body = np.vstack((sun_body_x, sun_body_y, sun_body_z))

    # Normal vector of panels in Body frame
    n_pan_body_x = np.einsum('ij,ij->j', n_pan, X_sat)
    n_pan_body_y = np.einsum('ij,ij->j', n_pan, Y_sat)
    n_pan_body_z = np.einsum('ij,ij->j', n_pan, Z_sat)
    n_pan_body = np.vstack((n_pan_body_x, n_pan_body_y, n_pan_body_z))

    # Modulation of solar flux including distance and eclipse
    FI_mod = ((astrody['UA'] / dr) ** 2 * nu)

    # Save results
    path_save = 'C:/Users/Carlo/OneDrive - INAF - Istituto Nazionale di Astrofisica/G4S/Matlab_ray/Results/'
    header_name = ['Satellite', 'T start [julian]', 'T End [julian]', 'Sampling time [day]']
    desc = pd.DataFrame({
        'VariableNames': header_name,
        'Values': [satellite, t_start, t_end, dt]
    })

    Ert_sat_dist = dr  # Earth-Sun distance

    # Cosine calculations
    cos_theta_sun = calculate_cos_theta_sun(r, X_sat, Y_sat, Z_sat, n_pan)

    # Force calculations
    Fx_sd, Fy_sd, Fz_sd = calculate_solar_forces(
        par_bw, cos_theta_sun, astrody, const, nu, dr, m_sat
    )

    # Satellite radiation pressure in Body frame
    SRP_body = np.vstack((Fx_sd, Fy_sd, Fz_sd)).T / m_sat
    time_utc = t

    # Save acceleration results
    np.savez(path_save + 'Sun_Body_ACC.npz', desc=desc, time_utc=time_utc, SRP_body=SRP_body)

    return init_Cond

# Helper Functions
def read_par(filename):
    params = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('%'):
                key_value = line.split('=')
                if len(key_value) == 2:
                    key = key_value[0].strip()
                    value = key_value[1].strip()
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                    params[key] = value
    return params

def read_box_wing_par(filename):
    # Implement parsing logic based on file structure
    par_bw = {}
    # ...
    return par_bw

def satellite_orb_ecc(a, om_t, MA_t, OMG_t, I, e):
    import numpy as np

    def solve_kepler(M, e):
        from scipy.optimize import newton
        func = lambda E: E - e * np.sin(E) - M
        func_prime = lambda E: 1 - e * np.cos(E)
        E_guess = M
        E = newton(func, E_guess, fprime=func_prime)
        return E

    E_t = np.array([solve_kepler(M, e) for M in MA_t])
    nu_t = 2 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E_t / 2),
        np.sqrt(1 - e) * np.cos(E_t / 2)
    )
    r = a * (1 - e * np.cos(E_t))

    x_orb = r * np.cos(nu_t)
    y_orb = r * np.sin(nu_t)

    cos_OMG = np.cos(OMG_t)
    sin_OMG = np.sin(OMG_t)
    cos_I = np.cos(I)
    sin_I = np.sin(I)
    cos_om = np.cos(om_t)
    sin_om = np.sin(om_t)

    r_sat = np.zeros((3, len(MA_t)))
    for i in range(len(MA_t)):
        R = np.array([
            [
                cos_OMG[i] * cos_om[i] - sin_OMG[i] * sin_om[i] * cos_I,
                -cos_OMG[i] * sin_om[i] - sin_OMG[i] * cos_om[i] * cos_I,
                sin_OMG[i] * sin_I
            ],
            [
                sin_OMG[i] * cos_om[i] + cos_OMG[i] * sin_om[i] * cos_I,
                -sin_OMG[i] * sin_om[i] + cos_OMG[i] * cos_om[i] * cos_I,
                -cos_OMG[i] * sin_I
            ],
            [
                sin_om[i] * sin_I,
                cos_om[i] * sin_I,
                cos_I
            ]
        ])
        r_orb = np.array([x_orb[i], y_orb[i], 0])
        r_sat[:, i] = R @ r_orb
    return r_sat

def umbra(r_sun, r_sat, RS, RT):
    # Implement the umbra function based on the geometric relations
    # between the Sun, Earth, and satellite positions
    nu = np.ones(r_sat.shape[1])  # Placeholder: no eclipse
    # ...
    return nu

def normalize_vectors(vectors):
    norms = np.linalg.norm(vectors, axis=0)
    return vectors / norms

def angle_between_vectors(v1, v2):
    v1_norm = normalize_vectors(v1)
    v2_norm = normalize_vectors(v2)
    sin_theta = np.linalg.norm(np.cross(v1_norm.T, v2_norm.T), axis=1)
    cos_theta = np.einsum('ij,ij->j', v1_norm, v2_norm)
    delta = np.arctan2(sin_theta, cos_theta)
    return delta

def modify_psi(Psi_nc, beta, eps, dt):
    rr = (np.abs(beta) < 4.1 / 180 * np.pi) & (np.abs(eps) < 10 / 180 * np.pi)
    rr_indices = np.where(rr)[0]
    n1 = np.diff(rr_indices)
    n2 = np.where(n1 != 1)[0]
    ind1 = np.concatenate(([0], n2 + 1))
    ind2 = np.concatenate((n2, [len(rr_indices) - 1]))
    fc = np.abs(eps[rr_indices[ind1 - 1]]) > 10 / 180 * np.pi
    ind = np.vstack((rr_indices[ind1[fc]], rr_indices[ind2[fc]]))

    Psi_meta = Psi_nc.copy()
    delta_Psi_meta = np.zeros_like(Psi_nc)
    for i in range(ind.shape[1]):
        idx_start = ind[0, i]
        idx_end = ind[1, i]
        old = Psi_nc[idx_start:idx_end + 1]
        sign_psi = np.sign(Psi_meta[idx_start - 1])
        Psi_meta[idx_start:idx_end + 1] = (
            90 * np.pi / 180 * sign_psi
            + (Psi_meta[idx_start - 1] - 90 * np.pi / 180 * sign_psi)
            * np.cos(2 * np.pi / 5656 * np.arange(1, idx_end - idx_start + 2) * dt * 86400)
        )
        delta_Psi_meta[idx_start:idx_end + 1] = old - Psi_meta[idx_start:idx_end + 1]
    return Psi_meta, delta_Psi_meta

def calculate_cos_theta_sun(r, X_sat, Y_sat, Z_sat, n_pan):
    cos_theta_sun = {}
    cos_theta_sun['X'] = np.einsum('ij,ij->j', r, X_sat)
    cos_theta_sun['Xp'] = np.clip(cos_theta_sun['X'], 0, None)
    cos_theta_sun['Xm'] = np.clip(-cos_theta_sun['X'], 0, None)
    cos_theta_sun['Y'] = np.einsum('ij,ij->j', r, Y_sat)
    cos_theta_sun['Yp'] = np.clip(cos_theta_sun['Y'], 0, None)
    cos_theta_sun['Ym'] = np.clip(-cos_theta_sun['Y'], 0, None)
    cos_theta_sun['Z'] = np.einsum('ij,ij->j', r, Z_sat)
    cos_theta_sun['Zp'] = np.clip(cos_theta_sun['Z'], 0, None)
    cos_theta_sun['Zm'] = np.clip(-cos_theta_sun['Z'], 0, None)
    cos_theta_sun['SA'] = np.einsum('ij,ij->j', r, n_pan)
    cos_theta_sun['SAp'] = np.clip(cos_theta_sun['SA'], 0, None)
    cos_theta_sun['SAm'] = np.clip(-cos_theta_sun['SA'], 0, None)
    cos_theta_sun['SA_X'] = np.einsum('ij,ij->j', X_sat, n_pan)
    cos_theta_sun['SAp_X'] = cos_theta_sun['SA_X']
    cos_theta_sun['SAm_X'] = -cos_theta_sun['SA_X']
    cos_theta_sun['SA_Y'] = np.einsum('ij,ij->j', Y_sat, n_pan)
    cos_theta_sun['SAp_Y'] = cos_theta_sun['SA_Y']
    cos_theta_sun['SAm_Y'] = -cos_theta_sun['SA_Y']
    cos_theta_sun['SA_Z'] = np.einsum('ij,ij->j', Z_sat, n_pan)
    cos_theta_sun['SAp_Z'] = cos_theta_sun['SA_Z']
    cos_theta_sun['SAm_Z'] = -cos_theta_sun['SA_Z']
    return cos_theta_sun

def calculate_solar_forces(par_bw, cos_theta_sun, astrody, const, nu, dr, m_sat):
    # Force calculations for different satellite faces
    Fxp, FXxp, FYxp, FZxp = calc_F_sun('Xp', par_bw, cos_theta_sun, astrody, const, nu, dr)
    Fxm, FXxm, FYxm, FZxm = calc_F_sun('Xm', par_bw, cos_theta_sun, astrody, const, nu, dr)
    Fyp, FXyp, FYyp, FZyp = calc_F_sun('Yp', par_bw, cos_theta_sun, astrody, const, nu, dr)
    Fym, FXym, FYym, FZym = calc_F_sun('Ym', par_bw, cos_theta_sun, astrody, const, nu, dr)
    Fzp, FXzp, FYzp, FZzp = calc_F_sun('Zp', par_bw, cos_theta_sun, astrody, const, nu, dr)
    Fzm, FXzm, FYzm, FZzm = calc_F_sun('Zm', par_bw, cos_theta_sun, astrody, const, nu, dr)
    FXsap, FYsap, FZsap = calc_F_SA_sun('SAp', par_bw, cos_theta_sun, astrody, const, nu, dr)
    FXsam, FYsam, FZsam = calc_F_SA_sun('SAm', par_bw, cos_theta_sun, astrody, const, nu, dr)

    Fx_sd = (
        Fxp - Fxm + 2 * (FXsap + FXsam) +
        FXxp + FXxm + FXyp + FXym + FXzp + FXzm
    )
    Fy_sd = (
        Fyp - Fym + 2 * (FYsap + FYsam) +
        FYxp + FYxm + FYyp + FYym + FYzp + FYzm
    )
    Fz_sd = (
        Fzp - Fzm + 2 * (FZsap + FZsam) +
        FZxp + FZxm + FZyp + FZym + FZzp + FZzm
    )
    return Fx_sd, Fy_sd, Fz_sd

def calc_F_sun(ax, par_bw, cos_theta, astrody, const, nu, d_sun_sat):
    F = 0
    FX = 0
    FY = 0
    FZ = 0
    Fp = 0
    for i in range(len(par_bw[ax])):
        com = (
            par_bw[ax][i][0] * (
                astrody['CS'] * (astrody['UA'] / d_sun_sat) ** 2 * nu
            ) * cos_theta[ax] / const['c']
        )
        F -= 2 * (
            (par_bw[ax][i][3]) / 3 + par_bw[ax][i][2] * cos_theta[ax]
        ) * com
        Fp -= com * (par_bw[ax][i][1] + par_bw[ax][i][3])
    FX = Fp * cos_theta['X']
    FY = Fp * cos_theta['Y']
    FZ = Fp * cos_theta['Z']
    return F, FX, FY, FZ

def calc_F_SA_sun(ax, par_bw, cos_theta, astrody, const, nu, d_sun_sat):
    F = 0
    Fp = 0
    for i in range(len(par_bw[ax])):
        com = (
            par_bw[ax][i][0] * (
                astrody['CS'] * (astrody['UA'] / d_sun_sat) ** 2 * nu
            ) * cos_theta[ax] / const['c']
        )
        F -= 2 * (
            (par_bw[ax][i][3]) / 3 + par_bw[ax][i][2] * cos_theta[ax]
        ) * com
        Fp -= com * (par_bw[ax][i][1] + par_bw[ax][i][3])
    FX = Fp * cos_theta['X'] + F * cos_theta[ax + '_X']
    FY = Fp * cos_theta['Y'] + F * cos_theta[ax + '_Y']
    FZ = Fp * cos_theta['Z'] + F * cos_theta[ax + '_Z']
    return FX, FY, FZ


# Questa funzione utilizza `m_sat`, `r`, `const`, `astrody`, e `par_bw` per calcolare SRP su ciascuna faccia
# del satellite in base alla distanza dal Sole e alla riflettività di ciascuna faccia.



def read_parameters(file_path):
    """
    Legge un file .par contenente costanti fisiche o astronomiche in formato chiave-valore.

    Args:
        file_path (str): Il percorso completo del file .par.

    Returns:
        dict: Dizionario con le costanti e i parametri letti dal file.
    """
    parameters = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Ignora righe vuote e commenti
            if not line or line.startswith('%'):
                continue
            # Divide la linea in chiave e valore
            if '=' in line:
                key, value = line.split('=')
                key = key.strip()
                value = value.split()[0].strip()  # Rimuove commenti aggiuntivi
                try:
                    # Converte il valore in float se possibile
                    parameters[key] = float(value)
                except ValueError:
                    parameters[key] = value  # Mantiene come stringa se non convertibile
    return parameters

def read_box_wing_parameters(file_path):
    """
    Legge i parametri del modello box-wing da un file .par, con facce come Xp, Xm, Yp, ecc.

    Args:
        file_path (str): Il percorso completo del file .par.

    Returns:
        dict: Dizionario con i parametri box-wing, per ciascuna faccia del satellite.
    """
    parameters = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Ignora righe vuote e commenti
            if not line or line.startswith('%'):
                continue
            # Divide la linea in chiave e valori
            if '=' in line:
                key, values_str = line.split('=')
                key = key.strip()
                values_str = values_str.split('%')[0].strip()  # Rimuove i commenti finali

                # Usa una regex per trovare i blocchi di numeri, gestendo sia "," che ";"
                try:
                    # Trova blocchi numerici separati da punto e virgola e convertili in liste di float
                    values = [
                        [float(num) for num in re.split(r'[,\s]+', block.strip())]
                        for block in values_str.strip('[];').split(';') if block
                    ]
                    parameters[key] = values
                except ValueError:
                    parameters[key] = values_str  # Mantiene come stringa se non convertibile
    return parameters


def mjuliandate(*args):
    # (Implementazione come sopra)


# %% main
# Definizione dei percorsi

 os.sys.path.append(r'C:/Users/Carlo/OneDrive - INAF - Istituto Nazionale di Astrofisica/G4S/Box_wing/Confronto Box Wing Surf Surf_Visco/Script di Confronto/Massimo_risultati')
# path_save = r'C:/Users/Carlo/OneDrive - INAF - Istituto Nazionale di Astrofisica/G4S/Matlab_ray/Results/'

path = r'C:/Users/carlo/OneDrive - INAF - Istituto Nazionale di Astrofisica/G4S/Matlab_ray'
sys.path.extend([path, os.path.join(path, 'Funzioni'), os.path.join(path, 'Dati Caomsol'), os.path.join(path, 'Class')])

Sat_inp = {"Box": "Foc3_Ansys_materials_9"}
wing = "yes"

Sat1 = SAT(Sat_inp["Box"], path)
Sat1_dat = read_mesh_mat_txt(Sat_inp["Box"], path)

if wing == "yes":
    Sat_inp["Wings"] = "Box_WING_Wings_center_009"

if len(Sat_inp) == 2:
    Sat2_data = read_mesh_mat_txt(Sat_inp["Wings"], path)
    Sat2 = SAT(Sat_inp["Wings"], path)

ray = c_Rays()

days = 0.5  
dminutes = 2  

deltaT = dminutes / (24 * 60)


# Ritestiamo la funzione completa
init_cond = Generazione_coordinate_Satellite_sole('E18', 0.01, 0.5)
init_cond
