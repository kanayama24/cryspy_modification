![cryspy_logo](https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/logo/cryspy_fix-03.png)

[![PyPI version](https://badge.fury.io/py/csp-cryspy.svg)](https://badge.fury.io/py/csp-cryspy)
[![Downloads](https://static.pepy.tech/badge/csp-cryspy)](https://pepy.tech/project/csp-cryspy)

# CrySPY
CrySPY (pronounced as crispy) is a crystal structure prediction tool written in Python.  
Document: https://tomoki-yamashita.github.io/CrySPY_doc  
Questions and comments: https://github.com/Tomoki-YAMASHITA/CrySPY/discussions

## Latest version
version 1.4.3 (2025 October 6)

## News
- [2025 Octorber 6] CrySPY 1.4.3 released.
    + Bug fix in EA-vc
- [2025 July 18] CrySPY 1.4.2 released.
    + Supports VASP and QE for EA-vc
    + Priority order of structure optimization input file names
    + Job file auto-rewriting
    + cryspy-skip subcommand
    + cryspy-calc-convex-hull subcommand
    + See [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info)
- [2025 July 7] CrySPY 1.4.1 released.
    + Support add_max, elim_max, and subs_max in EA-vc
    + Support charge neutral condition in EA-vc
    + New subcommand: cryspy-Eplot
    + See [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info)
- [2025 June 17] CrySPY 1.4.0 released.
    + Support variable-composition evolutionary algorithm
    + Support interactive mode with Jupyter
    + There are important changes. See [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info)
- [2024 May 31] CrySPY 1.3.0 released.
    + There are important changes. See [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info)
- [2024 May 10] CrySPY 1.2.5 released.
    + bug fix for order_ef in out_results.py
- [2024 May 7] CrySPY 1.2.4 released.
    + bug fix
- [2023 October 21] CrySPY 1.2.3 released.
    + bug fix for MPI
- [2023 October 18] CrySPY 1.2.2 released.
    + [Enthalpy](https://tomoki-yamashita.github.io/CrySPY_doc/features/enthalpy/index.html)
- [2023 September 27] CrySPY 1.2.1 released.
    + bug fix for ASE interface
- [2023 July 10] CrySPY 1.2.0 released. Version information/version 1.2.0
    + Interface for ASE
    + Adoption of logging
    + See also [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info) and [CHANGELOG](./CHANGELOG.md)
- [2023 June 14] CrySPY 1.1.1 released
    + bug fix
- [2023 May 16] CrySPY 1.1.0 released
    + MPI parallelization (optional)
    + New score of LAQA
    + See also [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info) and [CHANGELOG](./CHANGELOG.md)
- [2023 March 16] CrySPY 1.0.0 released
    + CrySPY is available in PyPI, so you can install by pip (project name is csp-cryspy).
    + See also [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info) and [CHANGELOG](./CHANGELOG.md)


## System requirements
### Python
- Python >= 3.8
- [PyXtal >= 0.5.3](https://pyxtal.readthedocs.io/en/latest "PyXtal")

(optional)
- [PHYSBO](https://www.pasums.issp.u-tokyo.ac.jp/physbo/en/about "PHYSBO") (required if algo is BO)
- [DScribe](https://singroup.github.io/dscribe/latest/ "DScribe") (required if algo is BO)
- [mpi4py](https://mpi4py.readthedocs.io/en/stable "mpi4py")
- [nglview](https://github.com/nglviewer/nglview "nglview")


See [CrySPY document](https://tomoki-yamashita.github.io/CrySPY_doc/installation/requirements/ "CrySPY document") in detail.

### Structure optimizer
At least one optimizer is required.

- [VASP](https://www.vasp.at "VASP") (tested with version 5.4.4)
- [QUANTUM ESPRESSO](http://www.quantum-espresso.org "Quantum ESPRESSO") (tested with version 6.x, version 5.x does not work)
- [OpenMX](http://www.openmx-square.org "OpenMX")
- [soiap](https://github.com/nbsato/soiap "soiap") (tested with version 0.2.2)
- [LAMMPS](http://lammps.sandia.gov "LAMMPS")
- [ASE](https://wiki.fysik.dtu.dk/ase "ASE")


## Document (English/Japanese)
[CrySPY document](https://tomoki-yamashita.github.io/CrySPY_doc "CrySPY documment")

## CrySPY Utility
[CrySPY Utility](https://github.com/Tomoki-YAMASHITA/CrySPY_utility "CrySPY Utility")

## Reference
### CrySPY (software)
* T. Yamashita, S. Kanehira, N. Sato, H. Kino, H. Sawahata, T. Sato, F. Utsuno, K. Tsuda, T. Miyake, and T. Oguchi, Sci. Technol. Adv. Mater. Meth. **1**, 87 (2021).
    - https://www.tandfonline.com/doi/full/10.1080/27660400.2021.1943171


### Bayesian optimization
* T. Yamashita, N. Sato, H. Kino, T. Miyake, K. Tsuda, and T. Oguchi, Phys. Rev. Mater. **2**, 013803 (2018).
    - https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.013803

* N. Sato, T. Yamashita, T. Oguchi, K. Hukushima, and T. Miyake, Phys. Rev. Mater. **4**, 033801 (2020).
    - https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.033801

### Baysian optimization and evolutionary algorithm
* T. Yamashita, H. Kino, K. Tsuda, T. Miyake, and T. Oguchi, Sci. Technol. Adv. Mater. Meth. **2**, 67 (2022).
    - https://www.tandfonline.com/doi/full/10.1080/27660400.2022.2055987

### LAQA
* K.Terayama, T. Yamashita, T. Oguchi, and K. Tsuda, npj Comput. Mater. **4**, 32 (2018).
    - https://www.nature.com/articles/s41524-018-0090-y

* T. Yamashita and H. Sekine, Sci. Technol. Adv. Mater. Meth. **2**, 84 (2022).
    - https://www.tandfonline.com/doi/full/10.1080/27660400.2022.2059335


## License
CrySPY is distributed under the MIT License.  
Copyright (c) 2018 CrySPY Development Team

#
# Cryspy - Slab Structure Exploration (by KK, 2025)

本バージョンのCryspyは、スラブ構造の探索に対応するよう改良されている。コード内の変更箇所には "modified by KK" または "added by KK" のコメントが付
されている。

## 拡張された入力パラメータ

### `struc_mode`
- 新たに `slab` モードを追加。
- `struc_mode=slab` を指定することで、スラブ構造探索モードが有効となる(selective dynamicsが`T T T`)。

### `r_relax`（float）
- `POSCAR_SLAB`における最表面層からバルク領域方向への距離（Å）を指定し、その範囲内の原子層を構造緩和の対象とする。

### `r_rearr`（float）
- `POSCAR_SLAB`における最表面層からバルク領域方向への距離（Å）を指定し、その範囲内の原子層を再構成の対象とする。
- 対象原子層は一度除去され、ランダムに原子が再配置される。
- 原子数は `nat` により指定（任意の原子数を指定可能）。

### `r_above`（float）
- `POSCAR_SLAB`における最表面層から真空領域方向への距離（Å）を指定し、ランダムに配置される原子のz座標の上限を示す。
- すなわち、再構成層の原子は最表面層から真空側に `r_above`、バルク側に `r_rearr` の範囲内でランダムに配置される。

### `dual_surface`（bool）
- `POSCAR_SLAB`の両面を再構成層とみなすか（`dual_surface=True`）、片面のみとするか（`dual_surface=False`）を指定する。
- `dual_surface=True` の場合、バルク部分の対称操作に基づき、片面の再構成層をもう一方に再配置する。
- 対称操作は以下のような行列で表される。ここで `*` は任意の値を示す：

```
[[*, *, *],
 [*, *, *],
 [0, 0, -1]]
```



## VASPによる構造探索の設定

`struc_mode=slab` を用いた構造探索をVASPで実行するには、以下の入力ファイルを `calc_in` ディレクトリに配置する必要がある：

- `INCAR_1`: 計算条件
- `KPOINTS`: k点サンプリング
- `POSCAR_SLAB`: スラブ構造
- `POTCAR`: 擬ポテンシャル
- `run.sh`: 実行スクリプト
