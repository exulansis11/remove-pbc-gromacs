#!/usr/bin/env python
"""Recupera os frames intactos de um XTC corrompido (inclusive com setores ilegiveis no disco).

Uso: xtc_scan.py ARQUIVO.xtc NATOMS [--extrai SAIDA.xtc] [--b T0] [--e T1]

O arquivo e lido em blocos; regioes que o disco nao consegue ler (erro de E/S) sao puladas.
Cada frame e localizado pelo cabecalho XTC (big-endian: magic 1995, natoms, step, time, box,
natoms) e o tamanho exato vem do campo de bytes comprimidos, entao so frames completos e
legiveis sao copiados. Sem --extrai, so lista os trechos continuos encontrados.
"""
import argparse
import os
import struct

p = argparse.ArgumentParser()
p.add_argument("xtc"); p.add_argument("natoms", type=int)
p.add_argument("--extrai"); p.add_argument("--b", type=float, default=-1); p.add_argument("--e", type=float, default=1e15)
a = p.parse_args()

SIG = struct.pack(">ii", 1995, a.natoms)
CAB = 92                     # cabecalho + parametros de compressao ate o n de bytes
BLOCO, PEDACO = 64 << 20, 1 << 20
fd = os.open(a.xtc, os.O_RDONLY)
tam = os.fstat(fd).st_size
ruins = 0                    # bytes ilegiveis


def le(pos, n):
    try:
        return os.pread(fd, n, pos)
    except OSError:
        return None


def frame_em(pos):
    """(tempo, tamanho) se ha um frame XTC valido e legivel em pos, senao None."""
    h = le(pos, CAB)
    if h is None or len(h) < CAB or h[:8] != SIG:
        return None
    t = struct.unpack(">f", h[12:16])[0]
    box = struct.unpack(">9f", h[16:52])
    if struct.unpack(">i", h[52:56])[0] != a.natoms or not 0.1 < box[0] < 1000 or t < 0:
        return None
    nb = struct.unpack(">i", h[88:92])[0]
    n = CAB + (nb + 3) // 4 * 4
    if nb <= 0 or pos + n > tam:
        return None
    return t, n


frames = []                  # (offset, tempo, tamanho)
pos = 0
while pos < tam:
    buf = le(pos, BLOCO + CAB)
    pedacos = [(pos, buf)] if buf is not None else []
    if buf is None:          # bloco com erro de leitura: tenta de 1 em 1 MB
        for q in range(pos, min(pos + BLOCO, tam), PEDACO):
            b = le(q, PEDACO + CAB)
            if b is None:
                ruins += PEDACO
            else:
                pedacos.append((q, b))
    for base, b in pedacos:
        i = b.find(SIG)
        while i != -1:
            off = base + i
            if not frames or off >= frames[-1][0] + frames[-1][2]:
                fr = frame_em(off)
                if fr:
                    frames.append((off, fr[0], fr[1]))
            i = b.find(SIG, i + 4)
    pos += BLOCO

print(f"{len(frames)} frames completos encontrados em {tam / 1e9:.2f} GB"
      + (f" | {ruins / 1e6:.0f} MB ilegiveis no disco (erro de E/S)" if ruins else ""))
if frames:
    dt = sorted(frames[k + 1][1] - frames[k][1] for k in range(min(len(frames) - 1, 50)))
    passo = dt[len(dt) // 2] if dt else 0
    ini = 0
    for k in range(1, len(frames) + 1):
        if k == len(frames) or frames[k][1] - frames[k - 1][1] > 1.5 * passo or frames[k][1] <= frames[k - 1][1]:
            print(f"  trecho: {frames[ini][1] / 1000:9.2f} -> {frames[k - 1][1] / 1000:9.2f} ns  ({k - ini} frames)")
            ini = k

if a.extrai:
    n = 0
    with open(a.extrai, "wb") as o:
        for off, t, nb in frames:
            if a.b - 1e-3 <= t <= a.e + 1e-3:
                d = le(off, nb)
                if d is not None and len(d) == nb:
                    o.write(d); n += 1
    print(f"extraidos {n} frames para {a.extrai}")
