#!/usr/bin/env python
"""Coloca cada cadeia da proteina e cada molecula de ligante na imagem periodica correta.

Uso: reimagem_ligante.py COMPLEXO.tpr TRAJ_SEM_FIT.xtc SAIDA.xtc [--contato 8]

Precisa de coordenadas ainda alinhadas com a caixa (antes de qualquer -fit). Por frame:
  1. Cadeias da proteina (fragmentos ligados): a partir da maior, cada cadeia vai para a imagem
     cujo centro fica mais perto do conjunto ja montado (mantem oligomeros inteiros).
  2. Cada molecula de ligante, separadamente: se alguma imagem encosta na proteina (< --contato A),
     usa a mais proxima; se esta solta no solvente, mantem a imagem continua com o frame anterior
     (assim um ligante solto nunca alterna de imagem) e so troca quando encosta numa copia
     periodica da proteina.
So translacoes exatas de vetores da caixa: a forma das moleculas nunca muda.
"""
import argparse
import itertools
import sys
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.lib.mdamath import triclinic_vectors

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("tpr"); ap.add_argument("xtc"); ap.add_argument("out")
ap.add_argument("--contato", type=float, default=8.0)
a = ap.parse_args()

u = mda.Universe(a.tpr, a.xtc)
prot, lig = u.select_atoms("protein"), u.select_atoms("not protein")
cadeias = sorted(prot.fragments, key=lambda f: -f.n_atoms)
ligantes = list(lig.fragments)
shifts = np.array(list(itertools.product(range(-2, 3), repeat=3)))

# Entrada ja alinhada (-fit)? Entao a caixa nao acompanha mais as coordenadas e transladar
# por vetores da caixa criaria contatos falsos: nesse caso, copia sem mexer.
ca = prot.select_atoms("name CA")
u.trajectory[0]
r0 = ca.positions - ca.center_of_geometry()
rot = 0.0
for ts in u.trajectory[::10]:
    R, _ = rotation_matrix(ca.positions - ca.center_of_geometry(), r0)
    rot = max(rot, np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
if rot < 3:
    with mda.Writer(a.out, u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(u.atoms)
    print(f"AVISO: entrada ja alinhada (rotacao max {rot:.1f} graus): reimagem nao aplicada")
    sys.exit(0)


def melhor_imagem(mol, alvo, vec, prev):
    """indice do vetor de vec para a molecula mol: contato com alvo, senao continuidade."""
    c, ca_ = mol.center_of_mass(), alvo.center_of_geometry()
    raio = np.linalg.norm(alvo.positions - ca_, axis=1).max() + np.linalg.norm(mol.positions - c, axis=1).max()
    dcom = np.linalg.norm(c + vec - ca_, axis=1)
    cand = np.where(dcom <= dcom.min() + 2 * raio)[0]
    cand = cand[np.argsort(dcom[cand])][:8]          # so as poucas imagens perto da proteina
    dmin = np.array([distance_array(alvo.positions, mol.positions + vec[k]).min() for k in cand])
    k = cand[dmin.argmin()]
    if dmin.min() > a.contato and prev is not None:  # solto: continuidade
        k = np.linalg.norm(c - ca_ + vec - prev, axis=1).argmin()
    return k, c - ca_ + vec[k]


prev = [None] * len(ligantes)
mov_cad = mov_lig = 0
with mda.Writer(a.out, u.atoms.n_atoms) as w:
    for ts in u.trajectory:
        vec = shifts @ triclinic_vectors(ts.dimensions)
        montado = cadeias[0]
        for cad in cadeias[1:]:
            d = np.linalg.norm(cad.center_of_geometry() + vec - montado.center_of_geometry(), axis=1)
            k = d.argmin()
            if np.linalg.norm(vec[k]) > 1e-3:
                cad.positions = cad.positions + vec[k]
                mov_cad += 1
            montado = montado + cad
        for i, mol in enumerate(ligantes):
            k, prev[i] = melhor_imagem(mol, prot, vec, prev[i])
            if np.linalg.norm(vec[k]) > 1e-3:
                mol.positions = mol.positions + vec[k]
                mov_lig += 1
        w.write(u.atoms)
print(f"{len(cadeias)} cadeia(s), {len(ligantes)} molecula(s) de ligante | "
      f"reposicionamentos: {mov_cad} de cadeia, {mov_lig} de ligante em {u.trajectory.n_frames} frames")
