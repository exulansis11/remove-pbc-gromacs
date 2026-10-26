#!/usr/bin/env bash
# remove_pbc.sh — remove o PBC de trajetorias proteina-ligante do GROMACS, verifica o resultado
# e (opcional) gera um PDB limpo frame a frame, sem PBC, sem saltos e sem ligante solto.
#
# Sem argumentos: abre um MENU que guia passo a passo.
#
# Modo direto:
#   remove_pbc.sh -s md.tpr -f md.xtc [-o prefixo] [-d pasta] [-b ps] [-e ps] [-dt ps]
#                 [-l RESNAME[,RESNAME...]] [--pdb] [--manter-soltos] [--min N]
#                 [--contato A] [--salto A]
#
#   -s  TPR da producao (sistema completo)       -f  XTC (bruto, "whole" ou corrompido)
#   -o  prefixo das saidas (padrao: nome do XTC)  -d  pasta de saida (padrao: atual)
#   -b/-e  janela em ps (padrao: trajetoria inteira)        -dt  passo em ps (padrao: 100)
#   -l  resname(s) do ligante (padrao: tudo que nao e proteina, agua ou ion)
#   --pdb            gera o XTC/PDB limpo (remove frames com PBC, saltos e ligante solto)
#   --manter-soltos  nao remove os frames com o ligante solto
#   --min N          trecho minimo de frames bons seguidos (padrao 10)
#   --contato A      ligante "solto" se a distancia minima a proteina for > A (padrao 6)
#   --salto A        "salto" se algum atomo andar > A entre frames seguidos (padrao 25)
#
# Autor: Genesio Martins de Aguiar Neto — https://github.com/exulansis11/remove-pbc-gromacs (MIT)
#
# Variaveis de ambiente: GMX (executavel do GROMACS), PYTHON (python com MDAnalysis),
#                        PBC_SOLVENTE (resnames de agua/ions; padrao abaixo).
set -euo pipefail
export LC_ALL=C
VERSAO=1.0.0
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LIB=$DIR/lib

# ------------------------------------------------------------------ configuracao
SOLV=${PBC_SOLVENTE:-"SOL WAT HOH TIP3 TIP3P TIP4 TIP4P TIP5P SPC SPCE T3P T4P NA CL K MG CA ZN NA+ CL- K+ CLA SOD POT"}
B=""; E=""; DT=100; PDB=nao; SOLTOS=remover; MIN=10; CONTATO=6; SALTO=25; LIGRES=""; PASTA=""; PFX=""

acha_gmx() {
  [ -n "${GMX:-}" ] && { command -v "$GMX" > /dev/null && return; echo "ERRO: GMX=$GMX nao encontrado"; exit 1; }
  for g in gmx gmx_mpi gmx_d gmx_mpi_d; do command -v $g > /dev/null && { GMX=$g; return; }; done
  echo "ERRO: GROMACS nao encontrado. Ative o ambiente (source GMXRC ou conda activate) ou defina GMX=..."
  exit 1
}
acha_python() {
  for p in ${PYTHON:-} python3 python; do
    [ -n "$p" ] && command -v "$p" > /dev/null && "$p" -c "import MDAnalysis, pandas" 2> /dev/null && { PY=$p; return; }
  done
  echo "ERRO: nenhum python com MDAnalysis e pandas. Instale: pip install -r requirements.txt"
  echo "      (ou defina PYTHON=/caminho/do/python)"
  exit 1
}
sel_solv() { local s; s=$(printf '"%s" ' $SOLV); echo "resname $s"; }   # resnames entre aspas (NA+ etc.)

linha() { printf '%*s\n' 62 '' | tr ' ' '='; }
log() { echo "[$PFX] $*"; }
q() { "$@" -quiet > "$TMP/gmx.log" 2>&1 || { echo "ERRO no GROMACS:"; tail -20 "$TMP/gmx.log"; exit 1; }; }
pergunta() {  # pergunta "texto" "padrao" -> RESP
  local p="$1" d="${2:-}"
  if [ -t 0 ]; then read -r -e -p "  $p${d:+ [$d]}: " RESP; else read -r RESP || RESP=""; fi
  RESP=${RESP:-$d}
}
sim_nao() { pergunta "$1 (s/n)" "$2"; [[ ${RESP,,} == s* || ${RESP,,} == y* ]]; }
arquivo() {   # arquivo "texto" -> RESP = caminho absoluto existente
  while true; do
    pergunta "$1"; RESP=${RESP/#\~/$HOME}
    [ -f "$RESP" ] && { RESP=$(realpath "$RESP"); return; }
    echo "  -> arquivo nao encontrado: $RESP"; [ -t 0 ] || exit 1
  done
}
entra_pasta() {
  [ -n "$1" ] || return 0
  mkdir -p "$1" && cd "$1" || { echo "ERRO: nao foi possivel usar a pasta $1"; exit 1; }
}
pyv() { "$PY" -c "import warnings; warnings.filterwarnings('ignore'); $1"; }

# lista os residuos nao-proteina e nao-solvente do TPR (candidatos a ligante)
candidatos() {
  pyv "
import MDAnalysis as m, collections
u = m.Universe('$1'); s = u.select_atoms('not protein and not resname ' + ' '.join('$SOLV'.split()))
c = collections.Counter(); r = collections.Counter()
for res in s.residues: c[res.resname] += res.atoms.n_atoms; r[res.resname] += 1
print(' '.join(f'{k}:{r[k]}:{c[k]}' for k in c))"
}
define_selecoes() {  # define LIGSEL (sintaxe do gmx select) a partir de LIGRES
  if [ -z "$LIGRES" ]; then   # automatico: todos os residuos que nao sao proteina nem solvente
    local c; c=$(candidatos "$TPR")
    [ -n "$c" ] || { echo "ERRO: nenhum residuo alem de proteina e solvente no TPR"; exit 1; }
    LIGRES=$(for x in $c; do echo -n "${x%%:*},"; done); LIGRES=${LIGRES%,}
  fi
  LIGSEL="resname $(echo "$LIGRES" | tr ',' ' ')"
}

# ------------------------------------------------------------------ processamento -> XTC
processa() {
  define_selecoes
  local jan=""; [ -n "$B$E" ] && jan="_${B:-0}-${E:-fim}ps"
  NOME=${PFX}_nopbc${jan}_dt${DT}
  local be=(); [ -n "$B" ] && be+=(-b "$B"); [ -n "$E" ] && be+=(-e "$E")

  log "1/5 conferindo integridade do XTC"
  # "file possibly corrupted" ou "Incomplete frame": em geral so um trecho esta danificado
  if grep -qiE "corrupt|incomplete frame" <<< "$($GMX check -f "$XTC" 2>&1 || true)"; then
    local nat; nat=$(pyv "import MDAnalysis as m; print(m.Universe('$TPR').atoms.n_atoms)")
    log "    XTC corrompido: recuperando os frames intactos"
    "$PY" "$LIB/xtc_scan.py" "$XTC" "$nat" --extrai "$TMP/rec.xtc" --b "${B:--1}" --e "${E:-1e15}" | sed "s/^/[$PFX]     /"
    XTC=$TMP/rec.xtc
  fi

  log "2/5 indices: proteina, ligante e complexo"
  q $GMX select -s "$TPR" -on "$TMP/sys.ndx" -select \
    "\"System\" all; \"Protein\" group \"Protein\"; \"Complex\" group \"Protein\" or ($LIGSEL); \"Ligand\" $LIGSEL"
  local nlig; nlig=$(awk '/\[ Ligand \]/{f=1;next} /\[/{f=0} f{n+=NF} END{print n+0}' "$TMP/sys.ndx")
  [ "$nlig" -gt 0 ] || { echo "ERRO: nenhum atomo de ligante com a selecao: $LIGSEL"; exit 1; }
  log "    ligante: $nlig atomos ($LIGSEL)"
  echo 2 | q $GMX convert-tpr -s "$TPR" -n "$TMP/sys.ndx" -o "${PFX}_complex.tpr"

  # entrada ja alinhada (-fit)? a caixa nao acompanha mais as coordenadas: nao refazer PBC
  local t0; t0=$($GMX dump -f "$XTC" 2> /dev/null | grep -m1 -o "time=[^ ]*" | cut -d= -f2 || true)
  t0=$(awk -v t="${t0:-0}" -v b="${B:-}" 'BEGIN{printf "%.0f", (b != "" ? b : t)}')
  echo 1 | q $GMX trjconv -s "$TPR" -f "$XTC" -n "$TMP/sys.ndx" -o "$TMP/amostra.pdb" -b "$t0" -e $((t0 + 3000)) -dt 100
  local rot; rot=$(pyv "
import numpy as np, MDAnalysis as m
from MDAnalysis.analysis.align import rotation_matrix
ca = m.Universe('$TMP/amostra.pdb').select_atoms('name CA'); ca.universe.trajectory[0]
r0 = ca.positions - ca.center_of_geometry(); a = 0
for ts in ca.universe.trajectory:
    R, _ = rotation_matrix(ca.positions - ca.center_of_geometry(), r0)
    a = max(a, np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
print(f'{a:.1f}' if ca.universe.trajectory.n_frames > 3 else '99')")

  if awk "BEGIN {exit !($rot < 3)}"; then
    log "3/5 AVISO: entrada ja alinhada (proteina girou ${rot} graus em 3 ns): PBC nao refeito"
    log "    para corrigir o PBC de verdade, use o XTC bruto da producao"
    echo 2 | q $GMX trjconv -s "$TPR" -f "$XTC" -n "$TMP/sys.ndx" -o "$TMP/recorte.xtc" "${be[@]}" -dt "$DT"
    log "4/5 desfazendo saltos do ligante solto"
    "$PY" "$LIB/desalta_ligante.py" "${PFX}_complex.tpr" "$TMP/recorte.xtc" "$TMP/reimg.xtc" | sed "s/^/[$PFX]     /"
  else
    log "3/5 whole -> nojump -> center/mol${B:+ (de $B ps)}${E:+ (ate $E ps)}, dt $DT ps"
    echo 0 | q $GMX trjconv -s "$TPR" -f "$XTC" -n "$TMP/sys.ndx" -o "$TMP/whole.xtc" -pbc whole "${be[@]}" -dt "$DT"
    echo 0 | q $GMX trjconv -s "$TPR" -f "$TMP/whole.xtc" -n "$TMP/sys.ndx" -o "$TMP/nojump.xtc" -pbc nojump
    printf '1\n2\n' | q $GMX trjconv -s "$TPR" -f "$TMP/nojump.xtc" -n "$TMP/sys.ndx" -o "$TMP/mol.xtc" -center -pbc mol
    log "4/5 ligante na imagem periodica correta"
    "$PY" "$LIB/reimagem_ligante.py" "${PFX}_complex.tpr" "$TMP/mol.xtc" "$TMP/reimg.xtc" | sed "s/^/[$PFX]     /"
  fi
  "$PY" "$LIB/diag_pbc.py" "$PFX" "${PFX}_complex.tpr" "$TMP/reimg.xtc" "${NOME}_diagnostico.csv" \
    --b -1 --e 1e15 --contato "$CONTATO" > /dev/null 2>&1
  log "5/5 alinhando pelo Backbone -> ${NOME}.xtc"
  q $GMX select -s "${PFX}_complex.tpr" -on "$TMP/cx.ndx" -select '"Complex" all; "Backbone" group "Backbone"'
  printf '1\n0\n' | q $GMX trjconv -s "${PFX}_complex.tpr" -f "$TMP/reimg.xtc" -n "$TMP/cx.ndx" -o "${NOME}.xtc" -fit rot+trans
}

# ------------------------------------------------------------------ verificacao
verifica() {  # verifica TPR TRAJ PREFIXO_CSV [DIAG]
  local tpr=$1 traj=$2 csv=$3 diag=${4:-}
  "$PY" "$LIB/audita_frames.py" "$PFX" "$tpr" "$traj" "${csv}_auditoria.csv" \
    --contato "$CONTATO" --salto "$SALTO" > /dev/null 2>&1
  "$PY" - "${csv}_auditoria.csv" "$diag" << 'EOF'
import sys, pandas as pd
a = pd.read_csv(sys.argv[1]); n = len(a); c = a.classe.value_counts()
pbc = 0
if sys.argv[2]:
    d = pd.read_csv(sys.argv[2]); pbc = int(d.classe.isin(["PBC_PROT", "PBC_LIG", "PBC_IMAGEM"]).sum())
q, s, so = int(c.get("QUEBRADO", 0)), int(c.get("SALTO", 0)), int(c.get("SOLTO", 0))
print(f"   frames analisados ............ {n}")
print(f"   com PBC (imagem/molecula) .... {pbc + q}")
print(f"   com salto entre frames ....... {s}")
print(f"   com ligante solto ............ {so}")
print(f"   ligacao mais longa ........... {a.lig_max_bond.max():.2f} A (ligante) / {a.prot_max_bond.max():.2f} A (proteina)")
if pbc + q == 0 and s == 0:
    v = "SEM PBC E SEM SALTOS" + (f" ({so} frames com ligante solto: afastamento fisico)" if so else "")
else:
    v = "ATENCAO: ha frames ruins — gere o XTC/PDB limpo para remove-los"
print(f"   VEREDITO: {v}")
EOF
}

# ------------------------------------------------------------------ XTC/PDB limpo
gera_limpo() {  # gera_limpo TPR_COMPLEXO XTC NOME [DIAG]
  local tpr=$1 xtc=$2 nome=$3 diag=${4:-}
  local opt=(); [ "$SOLTOS" = manter ] && opt+=(--manter-soltos); [ -n "$diag" ] && opt+=(--diag "$diag")
  log "auditoria frame a frame e remocao dos frames ruins"
  "$PY" "$LIB/audita_frames.py" "$PFX" "$tpr" "$xtc" "${nome}_auditoria.csv" \
    --contato "$CONTATO" --salto "$SALTO" > /dev/null 2>&1
  "$PY" "$LIB/limpa_frames.py" "$tpr" "$xtc" "${nome}_auditoria.csv" "${nome}_limpo" \
    --min "$MIN" --salto "$SALTO" "${opt[@]}" | sed "s/^/[$PFX]     /"
  log "gerando PDB (com CONECT) e verificando cada arquivo limpo"
  local tem_pymol=sim; "$PY" -c "import pymol" 2> /dev/null || tem_pymol=nao
  tail -n +2 "${nome}_limpo_blocos.csv" | while IFS=, read -r arq nf ini fim _; do
    echo 0 | q $GMX trjconv -s "$tpr" -f "$arq.xtc" -o "$arq.pdb" -conect
    echo "  ---- $arq.pdb ($nf frames, $ini-$fim ns)"
    verifica "$tpr" "$arq.xtc" "$arq"
    "$PY" "$LIB/checa_saltos.py" "$arq" "$arq.pdb" --limite "$SALTO" | sed 's/.*| proteina/   saltos no PDB: proteina/'
    [ $tem_pymol = sim ] && "$PY" "$LIB/checa_pymol.py" "$arq.pdb" "$tpr" 2> /dev/null | grep -v Warning | sed 's/^/   PyMOL: /'
  done
  return 0
}

resumo_final() {
  echo; linha; echo "  ARQUIVOS GERADOS (pasta: $(pwd))"; linha
  ls -1 "${PFX}"_* 2> /dev/null | sed 's/^/   /'
  echo
  if ls "${PFX}"_*_limpo.xtc > /dev/null 2>&1; then
    echo "  Use nas analises: *_limpo.xtc / *_limpo.pdb (sem PBC e sem saltos)."
    echo "  Frames retirados e motivo: *_limpo_frames_removidos.csv"
    ls "${PFX}"_*_limpo_bloco*.xtc > /dev/null 2>&1 && \
      echo "  *_limpo_blocoN: outros trechos bons, separados do principal por um salto."
  else
    echo "  XTC sem PBC: ${NOME:-}.xtc (para o PDB limpo, rode com --pdb ou use a opcao 3 do menu)."
  fi
  linha
}

# ------------------------------------------------------------------ menu
ajuda() {
  cat << 'TXT'

  COMO FUNCIONA
  1. Informe o TPR da producao e o XTC (bruto; se estiver corrompido, os frames intactos
     sao recuperados) e escolha qual residuo e o ligante.
  2. O PBC e removido: whole -> nojump -> center/mol -> ligante na imagem periodica correta
     -> alinhamento pelo Backbone. Sai um XTC so com proteina + ligante.
  3. O XTC e VERIFICADO: proteina e ligante inteiros, ligante na imagem certa, saltos.
  4. Se voce pedir o PDB, cada frame e auditado e os frames ruins sao REMOVIDOS (PBC,
     saltos e, se quiser, ligante solto). Os frames bons viram blocos sem salto, gravados
     em XTC + PDB com CONECT (sem ligacoes falsas no PyMOL), e tudo e verificado de novo.
     Nenhuma coordenada e alterada: so frames sao retirados.
  Detalhes: README.md

TXT
}

escolhe_ligante() {
  local c; c=$(candidatos "$TPR")
  if [ -z "$c" ]; then echo "  -> nenhum residuo alem de proteina e solvente no TPR"; exit 1; fi
  echo "  Residuos que nao sao proteina nem agua/ion:"
  for x in $c; do IFS=: read -r rn nr na <<< "$x"; printf '     %-8s %4s residuo(s), %6s atomos\n' "$rn" "$nr" "$na"; done
  local todos; todos=$(for x in $c; do echo -n "${x%%:*},"; done); todos=${todos%,}
  pergunta "Resname(s) do ligante, separados por virgula" "$todos"; LIGRES=$RESP
}

opcao_completa() {
  echo; echo "  >> REMOVER PBC (passo a passo)"; echo
  arquivo "TPR da producao (sistema completo)"; TPR=$RESP
  arquivo "XTC da producao"; XTC=$RESP
  escolhe_ligante
  local sug; sug=$(basename "$XTC" .xtc); sug=${sug%%[._]*}
  pergunta "Prefixo dos arquivos de saida" "$sug"; PFX=$RESP
  pergunta "Pasta de saida (sera criada)" "$PWD/${PFX}_sem_pbc"; PASTA=${RESP/#\~/$HOME}
  pergunta "Inicio da janela em ps (vazio = inicio da trajetoria)" "$B"; B=$RESP
  pergunta "Fim da janela em ps (vazio = fim da trajetoria)" "$E"; E=$RESP
  pergunta "Intervalo entre frames - dt (ps)" "$DT"; DT=$RESP
  echo; linha
  echo "  TPR ....... $TPR"; echo "  XTC ....... $XTC"; echo "  Ligante ... $LIGRES"
  echo "  Saida ..... $PASTA/${PFX}_*"; echo "  Janela .... ${B:-inicio} - ${E:-fim} ps, dt $DT ps"
  linha
  sim_nao "Confirmar e iniciar" s || return 0
  entra_pasta "$PASTA"
  roda_tudo pergunta_pdb
}

opcao_verificar() {
  echo; echo "  >> VERIFICAR TRAJETORIA OU PDB JA PROCESSADO"; echo
  arquivo "TPR do complexo (o *_complex.tpr gerado por este script)"; local tpr=$RESP
  arquivo "Trajetoria (.xtc ou .pdb)"; local tr=$RESP
  PFX=$(basename "${tr%.*}")
  entra_pasta "$(dirname "$tr")"
  echo; linha; echo "  VERIFICACAO de $(basename "$tr")"; linha
  verifica "$tpr" "$tr" "$PFX"
  "$PY" "$LIB/checa_saltos.py" "$PFX" "$tr" --limite "$SALTO" | sed 's/.*| proteina/   saltos: proteina/'
  if [[ $tr == *.pdb ]] && "$PY" -c "import pymol" 2> /dev/null; then
    "$PY" "$LIB/checa_pymol.py" "$tr" "$tpr" 2> /dev/null | grep -v Warning | sed 's/^/   PyMOL: /'
  fi
  linha; echo "  auditoria frame a frame: $(pwd)/${PFX}_auditoria.csv"
}

opcao_limpar() {
  echo; echo "  >> LIMPAR FRAMES RUINS DE UM XTC JA SEM PBC E GERAR O PDB"; echo
  arquivo "TPR do complexo (o *_complex.tpr gerado por este script)"; local tpr=$RESP
  arquivo "XTC sem PBC (so o complexo)"; local x=$RESP
  pergunta "Prefixo dos arquivos de saida" "$(basename "${x%.xtc}")"; PFX=$RESP
  pergunta "Pasta de saida (sera criada)" "$(dirname "$x")"; entra_pasta "${RESP/#\~/$HOME}"
  opcoes_limpeza
  gera_limpo "$tpr" "$x" "$PFX"
  resumo_final
}

opcoes_limpeza() {
  if sim_nao "Remover tambem os frames com o ligante SOLTO (> $CONTATO A da proteina)" s; then SOLTOS=remover; else SOLTOS=manter; fi
  pergunta "Trecho minimo de frames bons seguidos" "$MIN"; MIN=$RESP
}

roda_tudo() {   # $1 = pergunta_pdb (menu) ou vazio (modo direto)
  echo; linha; echo "  PROCESSANDO $PFX   (remove_pbc.sh $VERSAO)"; linha
  processa
  echo; linha; echo "  VERIFICACAO DO XTC (${NOME}.xtc)"; linha
  verifica "${PFX}_complex.tpr" "${NOME}.xtc" "$NOME" "${NOME}_diagnostico.csv"
  linha
  if [ "${1:-}" = pergunta_pdb ]; then
    echo
    if sim_nao "Gerar o PDB? (os frames ruins serao removidos antes)" s; then PDB=sim; opcoes_limpeza; fi
  fi
  [ "$PDB" = sim ] && gera_limpo "${PFX}_complex.tpr" "${NOME}.xtc" "$NOME" "${NOME}_diagnostico.csv"
  resumo_final
}

menu() {
  while true; do
    echo; linha
    echo "   REMOVE-PBC $VERSAO  —  proteina + ligante (GROMACS)"
    linha
    echo "   1) Remover PBC de uma trajetoria (passo a passo)"
    echo "   2) Verificar uma trajetoria ou PDB ja processado"
    echo "   3) Limpar frames ruins de um XTC ja sem PBC e gerar o PDB"
    echo "   4) Como funciona"
    echo "   0) Sair"
    linha
    pergunta "Escolha uma opcao" ""
    case $RESP in
      1) opcao_completa; return ;;
      2) opcao_verificar; return ;;
      3) opcao_limpar; return ;;
      4) ajuda ;;
      0 | q | Q | sair) return ;;
      *) echo "  opcao invalida"; [ -t 0 ] || return ;;
    esac
  done
}

# ------------------------------------------------------------------ entrada
uso() { sed -n '2,/^set -euo/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while [ $# -gt 0 ]; do
  case $1 in
    -s) TPR=$(realpath "$2"); shift 2 ;;  -f) XTC=$(realpath "$2"); shift 2 ;;
    -o) PFX=$2; shift 2 ;;  -d) PASTA=$2; shift 2 ;;  -l) LIGRES=$2; shift 2 ;;
    -b) B=$2; shift 2 ;;  -e) E=$2; shift 2 ;;  -dt) DT=$2; shift 2 ;;
    --pdb) PDB=sim; shift ;;  --manter-soltos) SOLTOS=manter; shift ;;
    --min) MIN=$2; shift 2 ;;  --contato) CONTATO=$2; shift 2 ;;  --salto) SALTO=$2; shift 2 ;;
    -v | --version) echo "remove_pbc.sh $VERSAO"; exit 0 ;;
    -h | --help) uso 0 ;;
    *) echo "opcao desconhecida: $1"; uso 1 ;;
  esac
done
acha_gmx; acha_python
TMP=$(mktemp -d "${TMPDIR:-/tmp}/remove_pbc.XXXXXX"); trap 'rm -rf "$TMP"' EXIT
if [ -z "${TPR:-}${XTC:-}" ]; then
  menu
  exit 0
fi
[ -n "${TPR:-}" ] && [ -n "${XTC:-}" ] || { echo "informe -s TPR e -f XTC"; uso 1; }
[ -n "$PFX" ] || { PFX=$(basename "$XTC" .xtc); PFX=${PFX%%[._]*}; }
entra_pasta "$PASTA"
roda_tudo
