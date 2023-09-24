# md_inputs

様々な分子の分子動力学シミュレーション(Molecular Dynamics; MD)向け入力ファイル

MDシミュレーションソフトウェアとしては、GENESIS、GROMACS, NAMD, Amberなどがありますが、ここではPythonに慣れた学生向けに主に[OpenMM](https://openmm.org)用の入力ファイルを置いています。

OpenMMは、minicondaをインストールした上で、以下のようにしたインストールしてください。

```
conda install -c conda-forge openmm
```

## 入力ファイル

- [adk_closed](adk_closed)
  -  Adenylate kinase の closed状態
- [adk_open](adk_open)
  - Adenylate kinase の open状態
- [lysozyme_R96H](lysozyme_R96H)
  - T4 lysozyme の変異体 (R96H)。熱安定性が天然より下がることがわかっている。
- [lysozyme_WT](lysozyme_WT)
  - T4 lysozyme の天然体
- [proteing](proteing)
  - Protein G に FRETの蛍光色素分子(donor & acceptor)をラベルしたもの
- [tau](tau)
  - Tau protein の PHF構造
