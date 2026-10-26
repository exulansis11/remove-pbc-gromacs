#!/usr/bin/env python
"""Confere um PDB multi-model do jeito que o PyMOL o le.

Uso: checa_pymol.py TRAJ.pdb [COMPLEXO.tpr]

Carrega o PDB no PyMOL (sem interface), pega as ligacoes que o PyMOL desenha e mede, em todos
os estados: a ligacao mais longa (> 2,5 A = parece PBC na tela) e o maior deslocamento de um
atomo entre estados consecutivos. Com o TPR, tambem compara as ligacoes do PyMOL com as reais.
Precisa do modulo pymol (ambiente md).
"""
import sys
import warnings

import numpy as np
from pymol import cmd

warnings.filterwarnings("ignore")
pdb = sys.argv[1]
cmd.load(pdb, "t")
n = cmd.count_states("t")
m = cmd.get_model("t", state=1)
b = np.array([x.index for x in m.bond])
ids = np.array([a.id - 1 for a in m.atom])
if len(sys.argv) > 2:
    import MDAnalysis as mda
    top = set(map(tuple, np.sort(mda.Universe(sys.argv[2]).bonds.indices, axis=1)))
    pym = set(map(tuple, np.sort(ids[b], axis=1)))
    print(f"ligacoes: PyMOL {len(pym)} | topologia {len(top)} | falsas {len(pym - top)} | faltando {len(top - pym)}")
longa, est_longa, salto, est_salto, prev = 0.0, 0, 0.0, 0, None
for st in range(1, n + 1):
    x = cmd.get_coords("t", state=st)
    L = np.linalg.norm(x[b[:, 0]] - x[b[:, 1]], axis=1).max()
    if L > longa:
        longa, est_longa = L, st
    if prev is not None:
        d = np.linalg.norm(x - prev, axis=1).max()
        if d > salto:
            salto, est_salto = d, st
    prev = x
print(f"{n} estados | ligacao mais longa {longa:.2f} A (estado {est_longa}) "
      f"{'OK' if longa < 2.5 else 'LIGACAO ESTICADA'} | maior deslocamento entre estados "
      f"{salto:.1f} A (estado {est_salto})")
