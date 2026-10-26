#!/usr/bin/env python
"""Procura saltos entre frames consecutivos numa trajetoria proteina-ligante (PDB multi-model ou TPR+XTC).

Uso: checa_saltos.py NOME TOPOLOGIA [TRAJ] [--csv saida.csv]

Por par de frames consecutivos (Angstrom):
  prot_max   maior deslocamento de um atomo da proteina
  lig_max    maior deslocamento de um atomo do ligante
  lig_com    deslocamento do centro de massa do ligante
Salto de PBC = um atomo atravessa a caixa (deslocamento da ordem do tamanho da caixa). Por
seguranca, qualquer deslocamento > --limite A (padrao 25) entre frames seguidos e marcado como
SALTO. Com dt = 100 ps, a difusao livre de um ligante pequeno fica em ~5-15 A, e atomos da
proteina alinhada raramente passam de ~15 A. Para dt maior, aumente --limite.
"""
import argparse
import warnings

import numpy as np
import MDAnalysis as mda

warnings.filterwarnings("ignore")
p = argparse.ArgumentParser()
p.add_argument("nome")
p.add_argument("top")
p.add_argument("traj", nargs="?")
p.add_argument("--csv")
p.add_argument("--limite", type=float, default=25.0)
a = p.parse_args()

u = mda.Universe(a.top, a.traj) if a.traj else mda.Universe(a.top)
prot = u.select_atoms("protein")
lig = u.select_atoms("not protein and not resname SOL HOH WAT NA CL K MG")
rows, prev = [], None
for ts in u.trajectory:
    cur = (prot.positions.copy(), lig.positions.copy(), lig.center_of_mass())
    if prev is not None:
        rows.append((ts.frame,
                     np.linalg.norm(cur[0] - prev[0], axis=1).max(),
                     np.linalg.norm(cur[1] - prev[1], axis=1).max(),
                     np.linalg.norm(cur[2] - prev[2])))
    prev = cur
r = np.array(rows)
if a.csv:
    np.savetxt(a.csv, r, delimiter=",", fmt=["%d", "%.3f", "%.3f", "%.3f"],
               header="frame,prot_max,lig_max,lig_com", comments="")
sp, sl = int((r[:, 1] > a.limite).sum()), int((r[:, 2] > a.limite).sum())
print(f"{a.nome}: {u.trajectory.n_frames} frames | proteina: desloc. max {r[:, 1].max():5.1f} A, "
      f"saltos>{a.limite:.0f}A = {sp} | ligante: desloc. max {r[:, 2].max():5.1f} A, "
      f"saltos>{a.limite:.0f}A = {sl} | {'SEM SALTOS' if sp + sl == 0 else 'COM SALTOS'}")
