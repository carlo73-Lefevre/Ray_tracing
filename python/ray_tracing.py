import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


@dataclass
class Ray:
    """Simple container for ray data."""
    p_source: Optional[np.ndarray] = None
    dir: Optional[np.ndarray] = None
    face_source: Optional[np.ndarray] = None


def facce_illuminate_potenziali(sat, ray) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """Identifica le facce potenzialmente illuminate.

    Parameters
    ----------
    sat : object
        Oggetto con gli attributi ``Normals_mesh`` e ``Centers_mesh``.
    ray : Ray
        Raggio sorgente con gli attributi ``p_source`` e ``dir``.

    Returns
    -------
    MAT : dict
        Dizionario con gli indici delle facce illuminate e restanti.
    loc : ndarray
        Indici delle facce illuminate nell'ordinamento per distanza.
    """
    # Prodotto scalare facce sorgenti
    ray_norm = ray.dir / np.linalg.norm(ray.dir, axis=1, keepdims=True)
    ddot = np.sum(sat.Normals_mesh * ray_norm, axis=1)
    idx_fc_ill = np.where(ddot < -1e-5)[0]

    # Calcola le distanze tra i centri delle facce e la sorgente di luce
    d = np.linalg.norm(ray.p_source - sat.Centers_mesh, axis=1)

    # Ordina le distanze e ottieni gli indici
    sorted_idx = np.argsort(d)

    # Trova le facce potenziali tra quelle rivolte nel verso giusto
    loc_mask = np.isin(sorted_idx, idx_fc_ill)

    # Filtra solo le facce potenzialmente illuminate
    facce_ill_idx = sorted_idx[loc_mask]
    facce_restanti = sorted_idx[~loc_mask]

    mat = {
        "Dis_sort_ill": facce_ill_idx,
        "facce_restanti": facce_restanti,
        "ordine_per_distanza": facce_ill_idx,
    }

    loc = np.nonzero(loc_mask)[0]
    return mat, loc


def dentro_list_vect_002(A: np.ndarray, B: np.ndarray, C: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Replica vettoriale della verifica barycentrica dei punti nei triangoli.

    Parameters
    ----------
    A, B, C : ndarray
        Vertici dei triangoli di riferimento (n x 3).
    P : ndarray
        Punti da verificare (n x 3). Si assume lo stesso numero di punti e triangoli.

    Returns
    -------
    ndarray
        Indici (0-based) dei triangoli che contengono un punto diverso dal proprio
        centro proiettato, ovvero le facce da escludere.
    """
    v0 = C - A
    v1 = B - A
    v2 = P[None, :, :] - A[:, None, :]

    dot00 = np.sum(v0 * v0, axis=1)[:, None]
    dot01 = np.sum(v0 * v1, axis=1)[:, None]
    dot11 = np.sum(v1 * v1, axis=1)[:, None]
    dot02 = np.sum(v0[:, None, :] * v2, axis=2)
    dot12 = np.sum(v1[:, None, :] * v2, axis=2)

    denom = dot00 * dot11 - dot01 * dot01
    u = (dot11 * dot02 - dot01 * dot12) / denom
    v = (dot00 * dot12 - dot01 * dot02) / denom

    n = A.shape[0]
    if n == P.shape[0]:  # escludi sovrapposizioni diagonali
        np.fill_diagonal(u, -1)
        np.fill_diagonal(v, -1)

    inside = (u > 0) & (v > 0) & (u + v < 1)
    tri_idx, pt_idx = np.where(inside)

    tri_ok = tri_idx.copy()
    mask = tri_idx < pt_idx
    tri_ok[mask] = pt_idx[mask]

    return np.unique(tri_ok)


def escludi_facce(sat, mat: Dict[str, np.ndarray], ray: Ray, cont: int) -> np.ndarray:
    """Esclude le facce in ombra utilizzando la proiezione sul piano sorgente."""

    idx_facce = mat["Dis_sort_ill"].copy()

    raggio_dir = ray.dir[cont]
    idx_p = sat.Faces[idx_facce, :]
    A = sat.Vertex[idx_p[:, 0], :]
    B = sat.Vertex[idx_p[:, 1], :]
    C = sat.Vertex[idx_p[:, 2], :]

    n_dir = np.tile(raggio_dir, (idx_facce.size, 1))

    t_f_a = np.sum(n_dir * (ray.p_source[cont] - A), axis=1) / np.sum(n_dir * n_dir, axis=1)
    t_f_b = np.sum(n_dir * (ray.p_source[cont] - B), axis=1) / np.sum(n_dir * n_dir, axis=1)
    t_f_c = np.sum(n_dir * (ray.p_source[cont] - C), axis=1) / np.sum(n_dir * n_dir, axis=1)

    F_A = A + t_f_a[:, None] * n_dir
    F_B = B + t_f_b[:, None] * n_dir
    F_C = C + t_f_c[:, None] * n_dir

    centri_mesh = sat.Centers_mesh[idx_facce, :]

    dentro_idx = dentro_list_vect_002(F_A, F_B, F_C, centri_mesh)
    if dentro_idx.size:
        idx_facce = np.delete(idx_facce, dentro_idx)

    return idx_facce


def ray_tracing(sat, sun_dir: np.ndarray, ray1_rif: np.ndarray, video: str = "off"):
    """Traduzione in Python di ``Ray_tracing.m``.

    Parameters
    ----------
    sat : object
        Oggetto con gli attributi ``Centers_mesh``, ``Faces``, ``Vertex`` e ``Normals_mesh``.
    sun_dir : ndarray
        Direzione del sole (3,).
    ray1_rif : ndarray
        Direzioni dei raggi riflessi per ciascuna faccia.
    video : str, optional
        Parametro ereditato dal codice MATLAB (non utilizzato).

    Returns
    -------
    sat : object
        Oggetto satellite (non modificato).
    rays : list of Ray
        Lista contenente il raggio sorgente e, se presenti, i raggi riflessi.
    frame : None
        Segnaposto per compatibilità con l'originale.
    """
    ray0 = Ray()

    n_faces = sat.Centers_mesh.shape[0]
    n_dir = np.tile(sun_dir, (n_faces, 1))
    cm = sat.Centers_mesh
    u_dir = np.tile(-sun_dir, (n_faces, 1))
    pn = np.tile(-sun_dir * 5, (n_faces, 1))

    t = np.sum(n_dir * (pn - cm), axis=1) / np.sum(n_dir * u_dir, axis=1)
    ray0.p_source = cm + t[:, None] * u_dir
    ray0.dir = -u_dir

    ray1 = Ray()

    mat, _ = facce_illuminate_potenziali(sat, ray0)
    ray0.p_source = ray0.p_source[mat["ordine_per_distanza"]]
    ray0.dir = ray0.dir[mat["ordine_per_distanza"]]

    mat["Dis_sort_ill"] = escludi_facce(sat, mat, ray0, 1)

    mat["source"] = mat["Dis_sort_ill"]
    mat["n_riflesso"] = np.ones(len(mat["Dis_sort_ill"]), dtype=int)

    all_faces = np.arange(sat.Faces.shape[0])
    mat["facce_scure"] = np.setdiff1d(all_faces, mat["facce_restanti"])

    if mat["source"].size > 0:
        ray1.p_source = sat.Centers_mesh[mat["source"], :]
        ray1.face_source = mat["source"]
        ray1.dir = ray1_rif[mat["source"], :]

    frame = None
    rays: List[Ray] = [ray0]
    if ray1.p_source is not None and ray1.p_source.size > 0:
        rays.append(ray1)

    return sat, rays, frame
