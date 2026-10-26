#!/usr/bin/env python
"""Auditoria frame a frame de uma trajetoria proteina-ligante ja sem PBC (XTC ou PDB).

Uso: audita_frames.py NOME COMPLEXO.tpr TRAJ.(xtc|pdb) SAIDA.csv [--contato 6] [--salto 25]
                      [--cadeias 8]

Por frame (Angstrom):
  lig_max_bond   ligacao mais longa do ligante (topologia)            > 2,5 = molecula cortada
  prot_max_bond  ligacao mais longa da proteina (topologia)           > 2,5 = molecula cortada
  cadeia_sep     maior distancia minima de uma cadeia as demais       > --cadeias = oligomero separado
  bb_rmsd        RMSD do backbone contra a estrutura media (informativo)
  mind           distancia minima a proteina da molecula de ligante MAIS AFASTADA; > --contato = SOLTO
  lig_desloc     maior deslocamento de um atomo do ligante vs frame anterior
  prot_desloc    maior deslocamento de um atomo da proteina vs frame anterior
Classe: QUEBRADO > SALTO (algum atomo anda > --salto) > SOLTO > OK
"""
import argparse
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.lib.distances import distance_array

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("nome"); ap.add_argument("tpr"); ap.add_argument("traj"); ap.add_argument("out")
ap.add_argument("--contato", type=float, default=6.0); ap.add_argument("--salto", type=float, default=25.0)
ap.add_argument("--cadeias", type=float, default=8.0)
arg = ap.parse_args()

u = mda.Universe(arg.tpr, arg.traj)
prot, lig = u.select_atoms("protein"), u.select_atoms("not protein")
cadeias = list(prot.fragments) if prot.n_atoms else []
ligantes = list(lig.fragments)
bb = u.select_atoms("backbone")
b = u.bonds.indices
pb, lb = b[np.isin(b[:, 0], prot.indices)], b[np.isin(b[:, 0], lig.indices)]
bbm = np.mean([bb.positions.copy() for ts in u.trajectory], axis=0)

rows, prev = [], None
for ts in u.trajectory:
    x = u.atoms.positions.copy()
    lmb = np.linalg.norm(x[lb[:, 0]] - x[lb[:, 1]], axis=1).max() if len(lb) else 0.0
    pmb = np.linalg.norm(x[pb[:, 0]] - x[pb[:, 1]], axis=1).max() if len(pb) else 0.0
    sep = 0.0
    for i, c in enumerate(cadeias if len(cadeias) > 1 else []):
        resto = sum((cadeias[j] for j in range(len(cadeias)) if j != i), cadeias[0][:0])
        sep = max(sep, distance_array(c.positions, resto.positions).min())
    r_bb = rmsd(bb.positions, bbm, center=True, superposition=True)
    mind = max(distance_array(prot.positions, m.positions).min() for m in ligantes)
    dl = np.linalg.norm(x[lig.indices] - prev[lig.indices], axis=1).max() if prev is not None else 0.0
    dp = np.linalg.norm(x[prot.indices] - prev[prot.indices], axis=1).max() if prev is not None else 0.0
    prev = x
    if lmb > 2.5 or pmb > 2.5 or sep > arg.cadeias:
        cls = "QUEBRADO"
    elif dl > arg.salto or dp > arg.salto:
        cls = "SALTO"
    elif mind > arg.contato:
        cls = "SOLTO"
    else:
        cls = "OK"
    rows.append((ts.frame, ts.time, lmb, pmb, sep, r_bb, mind, dl, dp, cls))

with open(arg.out, "w") as f:
    f.write("frame,tempo_ps,lig_max_bond,prot_max_bond,cadeia_sep,bb_rmsd,mind,lig_desloc,prot_desloc,classe\n")
    for r in rows:
        f.write("%d,%.0f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%s\n" % r)
cl = [r[-1] for r in rows]
print(f"{arg.nome}: frames={len(rows)} | OK={cl.count('OK')} SOLTO={cl.count('SOLTO')} "
      f"SALTO={cl.count('SALTO')} QUEBRADO={cl.count('QUEBRADO')} | {len(cadeias)} cadeia(s), "
      f"{len(ligantes)} molecula(s) de ligante")
