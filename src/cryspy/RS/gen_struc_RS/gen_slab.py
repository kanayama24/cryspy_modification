import numpy as np
from logging import getLogger
import random
from pymatgen.core import Structure
from pymatgen.analysis.local_env import VoronoiNN
from spglib import get_symmetry_dataset
#=========== modified by KK ============#
#from ...util.struc_util import check_distance, sort_by_atype
from ...util.struc_util import (
        check_distance, sort_by_atype,
        read_slab, symmetric_bottom, slab_site_properties
        )
#=======================================#
from pyxtal.tolerance import Tol_matrix

logger = getLogger('cryspy')

# coded by KK
def gen_slab_structure( 
        nstruc, 
        atype, 
        nat, 
        mindist,            
        r_relax, 
        r_rearr,
        r_above,
        dual_surface, 
        symprec=0.01,
        id_offset=0):

    '''
    Generate random structures for given slab structure

        tuple may be replaced by list

    # ---------- args
    nstruc (int): number of structures to be generated
    atype (tuple): atom type (e.g. ('Na', 'Cl'))
    nat (tuple): number of atoms (e.g. (4, 4)), None if vc=True
    mindist (): minimum interatomic distance (e.g. ((2.0, 1.5), (1.5, 2.0)))
    spgnum (str, int, or tuple): space group number 'all', 0, or tuple of space group numbers
    symprec (float): symmetry tolerance
    id_offset (int): structure ID starts from id_offset
                        e.g. nstruc = 3, id_offset = 10
                             you obtain ID 10, ID 11, ID 12

    # ---------- return
    init_struc_data (dict): {ID: pymatgen Structure, ...}
    '''

    # ---------- initialize
    init_struc_data = {}
    slab_init, slab_bulk, lattice, z_min, z_max  \
            = read_slab(r_rearr, dual_surface, atype)
    z_relax = r_relax/lattice.c
    z_rearr = r_rearr/lattice.c
    z_above = r_above/lattice.c

    # ---- structure generation
    while len(init_struc_data) < nstruc:
        tmp_struc = slab_bulk.copy()
        spg_sym_slab, spg_num_slab = tmp_struc.get_space_group_info(symprec=symprec)

        # ---- random generation of the top layer of the slab surface
        for elem, num in zip(atype, nat):
            for ___ in range(num):
                for __ in range(1000):
                    fpos = np.random.rand(3)
                    fpos[2] = np.random.uniform(z_max-z_rearr, z_max+z_above)
                    tmp_struc.append(elem, fpos, coords_are_cartesian=False)

                    # 距離のチェック
                    success, _, _ = check_distance(tmp_struc, atype, mindist)
                    if success:
                        break
                    else:
                        tmp_struc.pop()

        # ---- put atoms on bottom layer within symmety restriction
        if dual_surface:
            tmp_struc = symmetric_bottom(tmp_struc, slab_bulk, symprec=symprec)

        # ---- site properties (selective dynamics)
        tmp_struc = slab_site_properties(tmp_struc, dual_surface, z_max-z_relax, z_min+z_relax)

        # ---- check actual space group
        try:
            spg_sym, spg_num = tmp_struc.get_space_group_info(
                symprec=symprec)
        except TypeError:
            spg_num = 0
            spg_sym = None
        # ------ tmp_struc --> init_struc_data
        cid = len(init_struc_data) + id_offset
        init_struc_data[cid] = tmp_struc
        logger.info(f'Structure ID {cid:>6} was generated.'
                f' Space group: {spg_num_slab:>3} --> {spg_num:>3} {spg_sym}')

    return init_struc_data

