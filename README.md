# remove-pbc-gromacs

Remove as condições periódicas de contorno (PBC) de trajetórias **proteína–ligante** do GROMACS, **verifica** se o resultado ficou realmente sem PBC e, se você quiser, gera um **XTC/PDB limpo frame a frame**: sem PBC, sem saltos e sem ligante solto, pronto para SASA, superfície de contato, mapas de ligação de hidrogênio, MM-PBSA, visualização no PyMOL etc.

*English summary at the end.*

## Por que

O `gmx trjconv` sozinho (`-pbc mol -center`) costuma deixar problemas que estragam as análises sem dar sinal:

- o ligante fica na **imagem periódica errada**: longe da proteína nas coordenadas, mas encostado nela pela imagem mínima. SURF e H-bonds dão zero, e o MM-PBSA dá vdW = 0;
- com o ligante solto no solvente, a imagem **alterna a cada frame** e o ligante "salta" ~1 caixa;
- cadeias de um **oligômero** aparecem separadas;
- PDBs sem `CONECT` fazem o PyMOL **inventar ligações**, que esticam e parecem PBC;
- XTCs **corrompidos** param de ser lidos no meio.

Este script trata todos esses casos e mostra, frame a frame, o que sobrou.

## Instalação

Requer **GROMACS** (2018 ou mais novo; testado com 2021 e 2025) e **Python 3** com MDAnalysis e pandas.

```bash
git clone https://github.com/exulansis11/remove-pbc-gromacs.git
cd remove-pbc-gromacs
chmod +x remove_pbc.sh

# opção A: conda (instala GROMACS e Python juntos)
conda env create -f environment.yml
conda activate remove-pbc

# opção B: só o Python (se o GROMACS já está instalado)
pip install -r requirements.txt
```

A checagem no PyMOL é opcional: ela roda sozinha se o módulo `pymol` estiver instalado (`conda install -c conda-forge pymol-open-source`).

## Uso com menu

```bash
./remove_pbc.sh
```

```
==============================================================
   REMOVE-PBC 1.0.0  —  proteina + ligante (GROMACS)
==============================================================
   1) Remover PBC de uma trajetoria (passo a passo)
   2) Verificar uma trajetoria ou PDB ja processado
   3) Limpar frames ruins de um XTC ja sem PBC e gerar o PDB
   4) Como funciona
   0) Sair
==============================================================
```

Na **opção 1**, o script:
1. pede o TPR e o XTC da produção (a tecla Tab completa os caminhos);
2. lista os resíduos que não são proteína nem água ou íon, para você escolher qual é o ligante (útil quando há cofatores, como FAD ou heme);
3. pede o prefixo, a pasta de saída, a janela de tempo e o dt;
4. remove o PBC e **verifica o XTC**, mostrando um veredito (`SEM PBC E SEM SALTOS` ou `ATENCAO`);
5. **pergunta se é para gerar o PDB**. Se sim, audita frame a frame, **remove os frames ruins**, grava XTC + PDB limpos e verifica tudo de novo.

A **opção 2** só verifica um arquivo que já existe, sem alterar nada. A **opção 3** faz só a limpeza de frames, partindo de um XTC já sem PBC.

## Uso direto (scripts e lotes)

```bash
./remove_pbc.sh -s md.tpr -f md.xtc -o lig1 -d saida_lig1 -l LIG -b 100000 -e 170000 -dt 100 --pdb
```

| Opção | Significado | Padrão |
|---|---|---|
| `-s` | TPR da produção (sistema completo) | — |
| `-f` | XTC (bruto, "whole" ou corrompido) | — |
| `-o` | prefixo dos arquivos | nome do XTC |
| `-d` | pasta de saída | pasta atual |
| `-l` | resname(s) do ligante, separados por vírgula (ex.: `LIG` ou `GCG,FAD`) | tudo que não é proteína, água ou íon |
| `-b` / `-e` | janela em ps | trajetória inteira |
| `-dt` | intervalo entre frames em ps | 100 |
| `--pdb` | gera o XTC/PDB limpo | só XTC |
| `--manter-soltos` | não remove os frames com o ligante solto | remove |
| `--min N` | trecho mínimo de frames bons seguidos | 10 |
| `--contato A` | ligante "solto" se estiver a mais de A Å da proteína | 6 |
| `--salto A` | "salto" se algum átomo andar mais de A Å entre frames | 25 |

Variáveis de ambiente: `GMX` (executável do GROMACS, ex.: `gmx_mpi`), `PYTHON` (Python com MDAnalysis) e `PBC_SOLVENTE` (resnames de água e íons; veja o padrão no início do script).

Exemplo em lote:

```bash
for lig in lig1 lig2 lig3; do
  ./remove_pbc.sh -s $lig/md.tpr -f $lig/md.xtc -o $lig -d resultados/$lig -b 100000 -e 170000 --pdb
done
```

## Arquivos gerados

| Arquivo | Conteúdo |
|---|---|
| `X_nopbc_..._dt100.xtc` | trajetória sem PBC (só proteína + ligante), todos os frames da janela |
| `X_complex.tpr` | TPR só do complexo, para usar o XTC com `gmx` |
| `X_..._diagnostico.csv` | PBC frame a frame, antes do alinhamento |
| `X_..._auditoria.csv` | auditoria frame a frame do XTC |
| **`X_..._limpo.xtc` / `.pdb`** | **trajetória limpa, para as análises.** PDB com `CONECT` |
| `X_..._limpo_blocoN.xtc` / `.pdb` | outros trechos bons, separados do principal por um salto (ex.: o ligante mudou de sítio) |
| `X_..._limpo_frames_removidos.csv` | cada frame retirado e o motivo |
| `X_..._limpo_blocos.csv` | frames e tempo de cada arquivo limpo |

Na limpeza, **nenhuma coordenada é alterada**: só frames são retirados, e cada frame mantém o tempo original. Ao publicar, informe a fração de frames mantida (por exemplo, "587 de 701 frames; ligante ligado 84% do tempo").

## O que o script faz

| Etapa | Detalhe |
|---|---|
| 1. integridade | `gmx check`. Se o XTC estiver corrompido, `lib/xtc_scan.py` localiza os cabeçalhos dos frames intactos e os recupera (em geral, só um trecho do arquivo está danificado) |
| 2. índices | Protein, Ligand (pelo `-l` ou automático) e Complex, independentes de números de grupo |
| 3. PBC | `-pbc whole` → `-pbc nojump` → `-center` + `-pbc mol` |
| 4. imagens | `lib/reimagem_ligante.py`: remonta oligômeros (cada cadeia na imagem mais próxima das demais) e coloca **cada molécula** de ligante na imagem correta. Em contato com a proteína, usa a imagem que encosta; solta no solvente, mantém a continuidade com o frame anterior. Usa só translações exatas de vetores da caixa |
| 5. alinhamento | `-fit rot+trans` no Backbone, por último, porque depois dele a caixa não serve mais para corrigir imagens |
| 6. verificação | `lib/diag_pbc.py` (com a caixa, antes do fit) + `lib/audita_frames.py` (XTC final) |
| 7. limpeza (opcional) | `lib/limpa_frames.py` remove os frames ruins e separa os blocos sem salto; o PDB sai com `-conect`; cada arquivo é verificado de novo (`lib/checa_saltos.py`, `lib/checa_pymol.py`) |

**Entrada que já passou por `-fit`:** o script detecta isso sozinho (proteína com rotação < 3° em 3 ns) e não refaz o PBC, porque a caixa não acompanha mais as coordenadas e refazer pioraria o resultado. Ele só recorta e desfaz os saltos do ligante solto de forma aproximada (`lib/desalta_ligante.py`), sem nunca criar contatos. Para uma correção completa, use o XTC bruto.

## Frame ruim

| Motivo | Critério |
|---|---|
| `PBC` | proteína ou ligante cortados pela caixa, cadeias em imagens diferentes ou ligante na imagem periódica errada |
| `molecula quebrada` | ligação química > 2,5 Å ou cadeias do oligômero a > 8 Å umas das outras |
| `salto` | algum átomo se desloca mais que `--salto` Å entre frames seguidos |
| `ligante solto` | alguma molécula de ligante a mais de `--contato` Å da proteína (removido só se você escolher) |
| `trecho curto` | frames bons isolados entre frames ruins (menos que `--min` frames seguidos) |

"Ligante solto" não é erro de PBC: é o ligante se afastando de verdade, e isso também é um resultado (baixa afinidade, troca de sítio). Por isso o script pergunta se esses frames devem sair. Se o seu objetivo é estudar a dissociação, use `--manter-soltos`.

## Limitações

- Feito para **proteína + ligante(s)** em água. Membranas e ácidos nucleicos não foram testados.
- Proteínas sem contato entre si (por exemplo, duas proteínas livres na mesma caixa) seriam vistas como "oligômero separado".
- Os limites (`--contato`, `--salto`) foram calibrados para dt = 100 ps. Com dt bem maior, aumente `--salto`.
- O XTC sem PBC contém só o complexo (sem água). Análises que precisam do solvente, como pontes de água, devem usar a trajetória original.

## Testado em

- SARS-CoV-2 Mpro (monômero, 304 resíduos) com 6 ligantes, GROMOS 54A7, caixa dodecaédrica; inclui XTC corrompido, ligante que se solta e muda de sítio, e entradas já alinhadas.
- Tripanotiona redutase (dímero, 2 cadeias) com 2 FAD e 2 ligantes.

## Verificação automática (opcional)

O arquivo `extras/verifica.yml` confere, a cada envio, a sintaxe do script (`bash -n` e `shellcheck`), a compilação dos módulos Python e o `--help`. Para ativar no seu fork, copie-o para `.github/workflows/verifica.yml` (no GitHub: *Add file → Create new file*, cole o conteúdo e salve).

## Autor e licença

Genesio Martins de Aguiar Neto ([@exulansis11](https://github.com/exulansis11)). Distribuído sob a [licença MIT](LICENSE).

## Citação

Se usar, cite este repositório (https://github.com/exulansis11/remove-pbc-gromacs) e também o [GROMACS](https://www.gromacs.org) e o [MDAnalysis](https://www.mdanalysis.org) (Michaud-Agrawal et al., 2011; Gowers et al., 2016).

---

## English summary

`remove_pbc.sh` removes periodic boundary artifacts from GROMACS protein–ligand trajectories and **verifies** the result. It recovers corrupted XTC files, re-assembles oligomers, places every ligand molecule in the correct periodic image (nearest image when in contact, continuity when free in solvent), and fits on the backbone. Optionally (`--pdb`), it audits every frame and **removes bad frames** (PBC, jumps and, optionally, dissociated ligand), writing jump-free XTC/PDB blocks with `CONECT` records (so PyMOL does not invent bonds) and re-checking each output. No coordinates are modified during cleaning; only frames are dropped. Run it without arguments for an interactive menu (in Portuguese), or use the command-line options above. Requires GROMACS ≥ 2018 and Python 3 with MDAnalysis and pandas. MIT License, © Genesio Martins de Aguiar Neto.
