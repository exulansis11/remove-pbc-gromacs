#!/usr/bin/env python
"""Remove os frames ruins de uma trajetoria ja sem PBC e grava blocos sem saltos.

Uso: limpa_frames.py COMPLEXO.tpr TRAJ.xtc AUDITORIA.csv SAIDA [--diag DIAG.csv]
                     [--manter-soltos] [--min 10] [--salto 25]

Frame ruim:
  - QUEBRADO ou SALTO na auditoria (audita_frames.py)
  - PBC_PROT / PBC_LIG / PBC_IMAGEM no diagnostico de PBC (diag_pbc.py), se fornecido
  - SOLTO (ligante a > 6 A da proteina), a menos que --manter-soltos
Os frames bons formam trechos continuos; trechos com menos de --min frames sao descartados
(toques curtos). Trechos seguidos sao juntados num mesmo bloco se, entre o ultimo frame de um e o
primeiro do outro, nenhum atomo se desloca mais que --salto A; senao comecam um bloco novo
(ex.: ligante que mudou de sitio). Assim nenhum bloco tem salto.
Saidas: SAIDA.xtc (maior bloco), SAIDA_bloco2.xtc..., SAIDA_frames_removidos.csv, SAIDA_blocos.csv.
Nenhuma coordenada e alterada; o tempo original de cada frame e mantido.
"""
import argparse
import warnings

import numpy as np
import pandas as pd
import MDAnalysis as mda

warnings.filterwarnings("ignore")
p = argparse.ArgumentParser()
p.add_argument("tpr"); p.add_argument("xtc"); p.add_argument("auditoria"); p.add_argument("saida")
p.add_argument("--diag"); p.add_argument("--manter-soltos", action="store_true")
p.add_argument("--min", type=int, default=10); p.add_argument("--salto", type=float, default=25.0)
a = p.parse_args()

au = pd.read_csv(a.auditoria)
n = len(au)
motivo = np.array([""] * n, dtype=object)
for cl, txt in (("QUEBRADO", "molecula quebrada"), ("SALTO", "salto")):
    motivo[(au.classe == cl).values] = txt
if not a.manter_soltos:
    motivo[((au.classe == "SOLTO").values) & (motivo == "")] = "ligante solto (> 6 A)"
if a.diag:
    dg = pd.read_csv(a.diag)
    if len(dg) == n:
        pbc = dg.classe.isin(["PBC_PROT", "PBC_LIG", "PBC_IMAGEM"]).values
        motivo[pbc & (motivo == "")] = "PBC"
bom = motivo == ""

# trechos continuos de frames bons
trechos, i = [], 0
while i < n:
    if bom[i]:
        j = i
        while j + 1 < n and bom[j + 1]:
            j += 1
        if j - i + 1 >= a.min:
            trechos.append((i, j))
        else:
            motivo[i:j + 1] = f"trecho curto (< {a.min} frames)"
        i = j + 1
    else:
        i += 1

u = mda.Universe(a.tpr, a.xtc)
assert u.trajectory.n_frames == n, "auditoria e trajetoria com numeros de frames diferentes"


def pos(f):
    u.trajectory[f]
    return u.atoms.positions.copy()


# junta trechos seguidos sem salto na emenda
blocos = []
for t in trechos:
    if blocos and np.linalg.norm(pos(t[0]) - pos(blocos[-1][-1][1]), axis=1).max() <= a.salto:
        blocos[-1].append(t)
    else:
        blocos.append([t])
blocos.sort(key=lambda b: -sum(j - i + 1 for i, j in b))

linhas = []
for k, b in enumerate(blocos):
    nome = a.saida if k == 0 else f"{a.saida}_bloco{k + 1}"
    frames = [f for i, j in b for f in range(i, j + 1)]
    with mda.Writer(f"{nome}.xtc", u.atoms.n_atoms) as w:
        for f in frames:
            u.trajectory[f]
            w.write(u.atoms)
    linhas.append((nome, len(frames), au.tempo_ps[frames[0]] / 1000, au.tempo_ps[frames[-1]] / 1000, len(b)))
pd.DataFrame(linhas, columns=["arquivo", "frames", "inicio_ns", "fim_ns", "trechos"]).to_csv(
    f"{a.saida}_blocos.csv", index=False)
rem = pd.DataFrame({"frame": au.frame, "tempo_ns": au.tempo_ps / 1000, "motivo": motivo})
rem = rem[rem.motivo != ""]
rem.to_csv(f"{a.saida}_frames_removidos.csv", index=False)

print(f"frames: {n} | mantidos: {sum(l[1] for l in linhas)} | removidos: {len(rem)}")
for m, c in rem.motivo.value_counts().items():
    print(f"   removidos por {m}: {c}")
for nome, nf, ini, fim, nt in linhas:
    print(f"   {nome}.xtc: {nf} frames ({ini:.1f}-{fim:.1f} ns, {nt} trecho(s))")
