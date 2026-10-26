#!/usr/bin/env python
"""Diagnostico de PBC por frame para complexos proteina-ligante (coordenadas ANTES do fit).

Uso: diag_pbc.py NOME COMPLEXO.tpr TRAJ.xtc SAIDA.csv [--b ps] [--e ps] [--contato A]

Metricas por frame (Angstrom):
  prot_max_bond  ligacao mais longa da proteina (topologia)   > 2,5 = proteina cortada
  lig_max_bond   ligacao mais longa do ligante (topologia)    > 2,5 = ligante cortado
  mind_raw       distancia minima a proteina da molecula de ligante mais afastada (coordenadas)
  mind_mic       a mesma, pela imagem minima (distancia fisica real)
  cadeia_raw/mic maior distancia minima de uma cadeia as demais (coordenadas / imagem minima)
Classe:
  OK          tudo integro
  PBC_PROT    proteina cortada pela caixa, ou cadeias de um oligomero em imagens diferentes
  PBC_LIG     ligante cortado pela caixa
  PBC_IMAGEM  ligante longe nas coordenadas, mas encostado na proteina pela imagem minima
  DISSOCIADO  ligante longe mesmo pela imagem minima (afastamento fisico, nao PBC)
"""
import argparse
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

warnings.filterwarnings("ignore")
p = argparse.ArgumentParser()
p.add_argument("nome"); p.add_argument("top"); p.add_argument("traj"); p.add_argument("out")
p.add_argument("--b", type=float, default=-1); p.add_argument("--e", type=float, default=1e15)
p.add_argument("--contato", type=float, default=6.0)
a = p.parse_args()

u = mda.Universe(a.top, a.traj)
prot, lig = u.select_atoms("protein"), u.select_atoms("not protein")
if lig.n_atoms == 0:
    raise SystemExit(f"{a.nome}: ligante nao encontrado")
cadeias, ligantes = list(prot.fragments), list(lig.fragments)
b = u.bonds.indices
pb, lb = b[np.isin(b[:, 0], prot.indices)], b[np.isin(b[:, 0], lig.indices)]

def sep_cadeias(box):
    s = 0.0
    for i, c in enumerate(cadeias if len(cadeias) > 1 else []):
        resto = sum((cadeias[j] for j in range(len(cadeias)) if j != i), cadeias[0][:0])
        s = max(s, distance_array(c.positions, resto.positions, box=box).min())
    return s


rows = []
for ts in u.trajectory:
    if ts.time < a.b - 1e-3 or ts.time > a.e + 1e-3:
        continue
    x, box = u.atoms.positions, ts.dimensions
    has_box = box is not None and box[0] > 1
    pmb = np.linalg.norm(x[pb[:, 0]] - x[pb[:, 1]], axis=1).max() if len(pb) else 0.0
    lmb = np.linalg.norm(x[lb[:, 0]] - x[lb[:, 1]], axis=1).max() if len(lb) else 0.0
    mind_raw = max(distance_array(prot.positions, m.positions).min() for m in ligantes)
    mind_mic = (max(distance_array(prot.positions, m.positions, box=box).min() for m in ligantes)
                if has_box else np.nan)
    c_raw = sep_cadeias(None)
    c_mic = sep_cadeias(box) if has_box else np.nan
    if pmb > 2.5 or (has_box and c_raw > a.contato + 2 and c_mic <= a.contato + 2):
        cls = "PBC_PROT"
    elif lmb > 2.5:
        cls = "PBC_LIG"
    elif mind_raw > a.contato and has_box and mind_mic <= a.contato:
        cls = "PBC_IMAGEM"
    elif mind_raw > a.contato:
        cls = "DISSOCIADO"
    else:
        cls = "OK"
    rows.append((ts.time, pmb, lmb, mind_raw, mind_mic, c_raw, c_mic, cls))

with open(a.out, "w") as f:
    f.write("time_ps,prot_max_bond,lig_max_bond,mind_raw,mind_mic,cadeia_raw,cadeia_mic,classe\n")
    for r in rows:
        f.write("%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%s\n" % r)
cls = [r[-1] for r in rows]
print(f"{a.nome}: frames={len(rows)} | "
      + " ".join(f"{k}={cls.count(k)}" for k in ("OK", "PBC_PROT", "PBC_LIG", "PBC_IMAGEM", "DISSOCIADO")))
