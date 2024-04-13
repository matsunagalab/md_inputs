# Simulation of deca-alanine

- 1_build.py
  - simulation系のセットアップ。TerminalにACE、NMEでキャップを施し、C atom of the first peptide bondが原点に、 N atom of the last peptide bondがz軸に沿うように分子をorientする。

- 2_equilibration.py
  - Positional restraintsをCa atomsへ課して平衡化

- 3_production.py
  - Restrraint freeでの 100 ns simulation

- 3_umbrella.py
  - C atom of the first peptide bondを原点 (0, 0, 0)まわりに、 N atom of the last peptide bond を(0, 0, 1.4 nm)まわりでpositional restraintsをかけた 100 ns simulation


デフォルトではMac用のplathomeになっているので、Linuxで走らせる場合はLinuxのplathome部分のコメントを外してください(代わりにMac部分をコメント)。

