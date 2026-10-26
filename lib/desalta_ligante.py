#!/usr/bin/env python
"""Remove saltos de PBC das moleculas de ligante numa trajetoria que JA passou por -fit.

Uso: desalta_ligante.py COMPLEXO.tpr TRAJ_ALINHADA.xtc SAIDA.xtc [--contato 8]

Depois do fit, a caixa nao acompanha mais as coordenadas (a proteina girou antes do fit), entao
nao da para usar vetores da caixa fixos. Para cada molecula de ligante, cada salto (> 50 A em um
passo) e desfeito com o proprio vetor do salto, reescalado para o comprimento da caixa, e o
deslocamento acumulado e mantido nos frames seguintes (continuidade).
A molecula volta para a imagem original sempre que ela esta em contato com a proteina, e o
deslocamento nunca e aplicado se criaria um contato (o vetor e so estimado): os contatos nunca
mudam. Correcao aproximada: para corrigir de verdade, use o XTC bruto da producao.
"""
import argparse
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("tpr"); ap.add_argument("xtc"); ap.add_argument("out")
ap.add_argument("--contato", type=float, default=8.0)
a = ap.parse_args()
SALTO = 50.0

u = mda.Universe(a.tpr, a.xtc)
prot = u.select_atoms("protein")
ligantes = list(u.select_atoms("not protein").fragments)
desloc = [np.zeros(3) for _ in ligantes]
prev = [None] * len(ligantes)
movidos = saltos = 0
with mda.Writer(a.out, u.atoms.n_atoms) as w:
    for ts in u.trajectory:
        caixa = ts.dimensions[:3]
        for i, mol in enumerate(ligantes):
            rel = mol.center_of_mass() - prot.center_of_mass()
            if prev[i] is not None:
                d = rel + desloc[i] - prev[i]
                if np.linalg.norm(d) > SALTO:
                    L = caixa[np.abs(caixa - np.linalg.norm(d)).argmin()]
                    desloc[i] = desloc[i] - d / np.linalg.norm(d) * L
                    saltos += 1
            if np.linalg.norm(desloc[i]) > 1e-3:
                d_orig = distance_array(prot.positions, mol.positions).min()
                d_novo = distance_array(prot.positions, mol.positions + desloc[i]).min()
                if d_orig <= a.contato or d_novo <= a.contato:
                    desloc[i] = np.zeros(3)      # imagem original e a correta (contato)
                else:
                    mol.positions = mol.positions + desloc[i]
                    movidos += 1
            prev[i] = rel + desloc[i]
        w.write(u.atoms)
print(f"{len(ligantes)} molecula(s) de ligante | {saltos} saltos detectados; "
      f"{movidos} reposicionamentos (so com o ligante sem contato)")
