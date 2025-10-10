from collections import Counter
from logging import getLogger

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.periodic_table import DummySpecie

from ...util.struc_util import origin_shift, sort_by_atype, check_distance
from ...util.struc_util import (
        find_site, cal_g, sort_by_atype_mol, get_nat,
        read_slab, symmetric_bottom, slab_site_properties
        )

logger = getLogger('cryspy')


def gen_crossover_slab(
        atype,
        nat,
        mindist,
        struc_data,
        sp,
        n_crsov,
        id_start=None,
        symprec=0.01,
        nat_diff_tole=4,
        maxcnt_ea=50,
        r_relax=5,
        r_rearr=1,
        r_above=1,
        dual_surface=False,
    ):
    '''
    # ---------- args

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    nat (tuple): e.g. (4, 4)
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_crsov (int): number of structures to generate by crossover
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    nat_diff_tole (int): tolerance for nat_diff
    maxcnt_ea (int): maximum number of trial in crossover

    # ---------- return
    children (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, id of parent_B)}
    operation (dict): {id: 'crossover'}
    '''

    # ---------- initialize
    struc_cnt = 0
    children = {}
    parents = {}
    operation = {}

    # ---------- id_offset
    if id_start is None:
        cid = max(struc_data.keys()) + 1
    else:
        if id_start < (max(struc_data.keys()) + 1):
            logger.error('id_start is already included in structure ID of the data')
        else:
            cid = id_start

    # ---------- generate structures by crossover
    while struc_cnt < n_crsov:
        # ------ select parents
        pid_A, pid_B = sp.get_parents(n_parent=2)
        parent_A = struc_data[pid_A]
        parent_B = struc_data[pid_B]
        # ------ generate child
        child = gen_child(atype, nat, mindist, parent_A, parent_B, symprec,
                          nat_diff_tole, maxcnt_ea, r_relax, r_rearr, r_above, dual_surface)
        # ------ success
        if child is not None:
            children[cid] = child
            parents[cid] = (pid_A, pid_B)
            operation[cid] = 'crossover'
            try:
                spg_sym, spg_num = child.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            logger.info(f'Structure ID {cid:>6} was generated'
                    f' from {pid_A:>6} and {pid_B:>6} by crossover.'
                    f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return children, parents, operation


def gen_child(atype, nat, mindist, parent_A, parent_B, symprec,
              nat_diff_tole=4, maxcnt_ea=50, 
              r_relax=5, r_rearr=1, r_above=1, dual_surface=False):
    '''
    # ---------- args

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    nat (tuple): e.g. (4, 4)
    parent_A (Structure): pymatgen Structure object
    parent_B (Structure): pymatgen Structure object
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    nat_diff_tole (int): tolerance for nat_diff
    maxcnt_ea (int): maximum number of trial in crossover

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize
    #parent_A = origin_shift(parent_A)    # origin_shift returns a new Structure object
    #parent_B = origin_shift(parent_B)
    count = 0

    # ---------- separete surface layer from slab
    slab_init, slab_bulk, lattice, z_min, z_max \
            = read_slab(r_rearr, dual_surface, atype)
    z_mean = (z_min+z_max)/2
    nsite_bulk = len(slab_bulk.sites)
    bulk_species = [site.species_string for site in slab_bulk.sites]
    bulk_coords = [site.frac_coords for site in slab_bulk.sites]
    #======== to be deleted ===========#
    #bulk_coords = [(site_1.frac_coords+site_2.frac_coords)/2 for site_1, site_2 in zip(bulk_pA_sites, bulk_pB_sites)]
    #
    # ---------- check identity of bulk structure
    #bulk_pA_sites = parent_A.sites[:nsite_bulk]
    #bulk_pB_sites = parent_B.sites[:nsite_bulk]
    #if ([site.species_string for site in bulk_pA_sites] != bulk_species or 
    #    [site.species_string for site in bulk_pB_sites] != bulk_species):
    #    logger.error('parent_A or parent_B seem to have different bulk layers')
    #==================================#

    # ---------- surface layer structure
    parent_A_rearr_sites = parent_A.sites[nsite_bulk:]
    parent_B_rearr_sites = parent_B.sites[nsite_bulk:]
    parent_A_rearr = Structure(
        lattice=lattice,
        species=[site.species_string for site in parent_A_rearr_sites if site.frac_coords[2]>z_mean],
        coords=[site.frac_coords for site in parent_A_rearr_sites if site.frac_coords[2]>z_mean],
        coords_are_cartesian=False 
    )
    parent_B_rearr = Structure(
        lattice=lattice,
        species=[site.species_string for site in parent_B_rearr_sites if site.frac_coords[2]>z_mean],
        coords=[site.frac_coords for site in parent_B_rearr_sites if site.frac_coords[2]>z_mean],
        coords_are_cartesian=False
    )

    # ---------- generate child
    while True:
        count += 1
        # ------ coordinate crossover
        axis, slice_point, species, coords = _one_point_crossover(parent_A_rearr, parent_B_rearr)
        # ------ child structure
        child_rearr = Structure(lattice, species, coords)
        child = Structure(lattice, bulk_species+species, bulk_coords+coords)
        # ------ check nat_diff
        nat_diff = _get_nat_diff(atype, nat, child_rearr)
        if any([abs(n) > nat_diff_tole for n in nat_diff]):
            logger.debug(f'nat_diff = {nat_diff}')
            if count > maxcnt_ea:    # fail
                return None
            continue    # slice again
        # ------ check mindist
        success, _, _ = check_distance(child, atype, mindist, check_all=False)
        # ------ something smaller than mindist
        if not success:
            # -- remove atoms within mindist
            if any([n > 0 for n in nat_diff]):
                child = _remove_within_mindist(child, atype, mindist, nat_diff, nsite_bulk)
                if child is None:    # fail --> slice again
                    if count > maxcnt_ea:
                        return None
                    continue
            else:    # nothing to remove, nat_diff = [0, 0]
                if count > maxcnt_ea:
                    return None
                continue    # fail --> slice again
        # ------ recheck nat_diff
        # ------ excess of atoms
        nat_diff = _get_nat_diff(atype, nat, child_rearr)    # recheck
        if any([n > 0 for n in nat_diff]):
            child = _remove_border_line(child, atype, axis,
                                        slice_point, nat_diff, nsite_bulk)
        # ------ lack of atoms
        nat_diff = _get_nat_diff(atype, nat, child_rearr)    # recheck
        if any([n < 0 for n in nat_diff]):
            child = _add_border_line(child, atype, mindist, axis, slice_point,
                                     nat_diff, nsite_bulk, maxcnt_ea)
        # ------ success --> break while loop
        if child is not None:
            break
        # ------ fail --> slice again
        else:
            if count > maxcnt_ea:
                return None
            continue

    # ---------- final check for nat
    nat_diff = _get_nat_diff(atype, nat, child_rearr)
    if not all([n == 0 for n in nat_diff]):
        return None    # failure

    # ---------- put atoms on bottom layer within symmetry
    if dual_surface:
        child = symmetric_bottom(child, slab_bulk,symprec=symprec)

    # ---------- site properties (selective dynamics)
    z_relax = r_relax/lattice.c
    child = slab_site_properties(child, z_max-z_relax, z_min+z_relax)

    # ---------- return
    return child


def _one_point_crossover(parent_A, parent_B):
    # ---------- slice point
    while True:
        slice_point = np.random.normal(loc=0.5, scale=0.1)
        if 0.3 <= slice_point <= 0.7:
            break
    axis = np.random.choice([0, 1])

    # ---------- crossover
    species_A = []
    species_B = []
    coords_A = []
    coords_B = []
    for i in range(parent_A.num_sites):
        if parent_A.frac_coords[i, axis] <= slice_point:
            species_A.append(parent_A[i].species_string)
            coords_A.append(parent_A[i].frac_coords)
        else:
            species_B.append(parent_A[i].species_string)
            coords_B.append(parent_A[i].frac_coords)
    for i in range(parent_B.num_sites):
        if parent_B.frac_coords[i, axis] >= slice_point:
            species_A.append(parent_B[i].species_string)
            coords_A.append(parent_B[i].frac_coords)
        else:
            species_B.append(parent_B[i].species_string)
            coords_B.append(parent_B[i].frac_coords)

    # ---------- adopt a structure with more atoms
    if len(species_A) > len(species_B):
        species = species_A
        coords = coords_A
    elif len(species_A) < len(species_B):
        species = species_B
        coords = coords_B
    else:
        if np.random.choice([0, 1]):
            species = species_A
            coords = coords_A
        else:
            species = species_B
            coords = coords_B

    # ---------- return
    return axis, slice_point, species, coords


def _get_nat_diff(atype, nat, child):
    '''
    original nat - child nat
    e.g.
        nat = [4, 4]        # original
        tmp_nat = [3, 5]    # child
        nat_diff = [-1, 1]
    '''
    tmp_nat = get_nat(child, atype)
    nat_diff = [i - j for i, j in zip(tmp_nat, nat)]
    return nat_diff


def _remove_within_mindist(child, atype, mindist, nat_diff, nsite_bulk):
    '''
    if success: return child
    if fail:    return None
    '''
    for itype in range(len(atype)):
        while nat_diff[itype] > 0:
            # ---------- check dist
            dist_list = check_distance(child, atype, mindist, check_all=True)
            if not dist_list:    # nothing within mindist
                return child
            # ---------- appearance frequency
            ij_within_dist = [isite[0] for isite in dist_list if isite[0]>=nsite_bulk] + [
                jsite[1] for jsite in dist_list if jsite[1]>=nsite_bulk]
            site_counter = Counter(ij_within_dist)
            # ---------- get index for removing
            rm_index = None
            # ---- site[0]: index, site[1]: count
            for site in site_counter.most_common():
                if child[site[0]].species_string == atype[itype]:
                    rm_index = site[0]
                    break    # break for loop
            # ---------- remove atom
            if rm_index is None:    # fail
                return None
            else:
                child.remove_sites([rm_index])
                nat_diff[itype] -= 1

    # ---------- final check
    dist_list = check_distance(child, atype, mindist, check_all=True)
    if dist_list:    # still something within mindist
        logger.warning('remove_within_mindist: some atoms within mindist. retry.')
        return None
    else:    # success
        return child


def _remove_border_line(child, atype, axis, slice_point, nat_diff, nsite_bulk):
    # ---------- rank atoms from border line
    coords_axis = child.frac_coords[:, axis]

    # ---------- boundary --> 0.0, slice_point, 1.0
    near_sp = (slice_point/2.0 < coords_axis) & \
        (coords_axis < (slice_point + 1.0)/2.0)
    near_one = (slice_point + 1.0)/2.0 <= coords_axis

    # ---------- distance from nearest boundary
    coords_diff = np.where(near_sp,
                            abs(coords_axis - slice_point),
                            coords_axis)
    coords_diff = np.where(near_one, 1.0 - coords_diff, coords_diff)
    atom_border_indx = np.argsort(coords_diff)

    # ---------- remove list
    rm_list = []
    for itype, nrm in enumerate(nat_diff):
        rm_list.append([])
        if nrm > 0:
            for ab_indx in atom_border_indx:
                if child[ab_indx].species_string == atype[itype]:
                    if ab_indx >= nsite_bulk:
                        rm_list[itype].append(ab_indx)
                if len(rm_list[itype]) == nrm:
                    break

    # ---------- remove
    for each_type in rm_list:
        if each_type:
            child.remove_sites(each_type)

    # ---------- return
    return child


def _add_border_line(child, atype, mindist, axis, slice_point, nat_diff, nsite_bulk, maxcnt_ea=50):
    for i in range(len(atype)):
        # ---------- counter
        cnt = 0

        # ---------- add atoms
        while nat_diff[i] < 0:
            cnt += 1
            coords = np.random.rand(3)
            mean = _mean_choice(child, axis, slice_point)
            coords[axis] = np.random.normal(loc=mean, scale=0.08)
            coords[2] = np.random.normal(loc=np.mean([site.frac_coords[2] for site in child.sites[nsite_bulk:]]), scale=0.08)
            child.append(species=atype[i], coords=coords)
            success, mindist_ij, dist = check_distance(child, atype, mindist)
            if success:
                cnt = 0    # reset
                nat_diff[i] += 1
                print("add an atom")
            else:
                type0 = atype[mindist_ij[0]]
                type1 = atype[mindist_ij[1]]
                logger.warning(f'mindist in _add_border_line: {type0} - {type1}, {dist}. retry.')
                child.pop()    # cancel
            # ------ fail
            if cnt == maxcnt_ea:
                return None

        # ---------- return
        return child


def _mean_choice(child, axis, slice_point):
    '''
    Which border contains the most atoms?
    '''
    n_zero = np.sum(np.abs(child.frac_coords[:, axis] - 0.0)
                    < 0.1)
    n_slice = np.sum(np.abs(child.frac_coords[:, axis]
                            - slice_point) < 0.1)
    if n_zero < n_slice:
        mean = 0.0
    elif n_zero > n_slice:
        mean = slice_point
    else:
        mean = np.random.choice([0.0, slice_point])
    return mean
