import numpy as np
from logging import getLogger
import random
from pymatgen.core import Structure
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from spglib import get_symmetry_dataset
from ...util.struc_util import check_distance, sort_by_atype
from pyxtal.tolerance import Tol_matrix

logger = getLogger('cryspy')

# coded by KK
def gen_slab_structure( 
        nstruc, 
        atype, 
        nat, 
        mindist,            
        z_top, 
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
    tmp_nat = nat
    tmp_atype = atype
    tolmat = _set_tol_mat(tmp_atype, mindist)

    slab = Structure.from_file('./calc_in/POSCAR_SLAB')
    lattice = slab.lattice
    z_frac_max = max([site.frac_coords[2] for site in slab])
    z_frac_min = min([site.frac_coords[2] for site in slab])
    z_cut = z_frac_max + z_top / lattice.c

    def random_frac_in_top():  #random atom generation
        f = np.random.rand(3)
        f[2] = np.random.uniform(z_frac_max, z_cut)
        return f

    def decorate_top(tmp_struc): #random generation of the top layer of the slab surface
        for elem, num in zip(tmp_atype, tmp_nat):
            for ___ in range(num):
                for __ in range(1000):
                    fpos = random_frac_in_top()

                    # 既存の selective_dynamics を取得（なければ空リスト）
                    current_flags = tmp_struc.site_properties.get("selective_dynamics", [])

                    # 更新
                    tmp_struc.append(elem, fpos, coords_are_cartesian=False)
                    current_flags.append(np.array([True,True,True]))
                    tmp_struc.add_site_property("selective_dynamics", current_flags)

                    # 距離のチェック
                    success, _, _ = check_distance(tmp_struc, tmp_atype, mindist)
                    if success:
                        break
                    else:
                        tmp_struc.pop()
        return tmp_struc

    def symmetric_bottom(tmp_struc): #bottm layer generation by symmetry operation
        
        # 対称性情報の取得
        analyzer = SpacegroupAnalyzer(slab, symprec=symprec)
        dataset = analyzer.get_symmetry_dataset()
        for rot, trans in zip(dataset.rotations, dataset.translations):
            if np.allclose(rot[2], [0,0,-1]):
                break
        #print("symmop(trans):",trans)
        #print("symmop(rot):")
        #print(rot)

        # 対照操作による新しい原子の追加
        for site in tmp_struc.sites:
            fpos = site.frac_coords
            if fpos[2] > z_frac_max:
                new_fpos = np.dot(rot,fpos) + trans
                new_fpos = new_fpos % 1.0

                # 既存の selective_dynamics を取得（なければ空リスト）
                current_flags = tmp_struc.site_properties.get("selective_dynamics", [])

                # 更新
                tmp_struc.append(site.specie, new_fpos, coords_are_cartesian=False)
                current_flags.append(np.array([True,True,True]))
                tmp_struc.add_site_property("selective_dynamics", current_flags)
        return tmp_struc

    while len(init_struc_data) < nstruc:
        tmp_struc = slab.copy()
        spg_sym_slab, spg_num_slab = tmp_struc.get_space_group_info(
                                    symprec=symprec)
        tmp_struc = decorate_top(tmp_struc)
        if dual_surface == True:
            tmp_struc = symmetric_bottom(tmp_struc)
        else:
            pass

        # -- check actual space group
        try:
            spg_sym, spg_num = tmp_struc.get_space_group_info(
                symprec=symprec)
        except TypeError:
            spg_num = 0
            spg_sym = None
        # ------ tmp_struc --> init_struc_data
        #tmp_struc = sort_by_atype(tmp_struc, atype)
        cid = len(init_struc_data) + id_offset
        init_struc_data[cid] = tmp_struc
        logger.info(f'Structure ID {cid:>6} was generated.'
                f' Space group: {spg_num_slab:>3} --> {spg_num:>3} {spg_sym}')

    return init_struc_data

def _set_tol_mat(atype, mindist):
    tolmat = Tol_matrix()
    for i, itype in enumerate(atype):
        for j, jtype in enumerate(atype):
            if i <= j:
                tolmat.set_tol(itype, jtype, mindist[i][j])
    return tolmat

