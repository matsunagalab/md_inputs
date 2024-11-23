# README.md
実行手順は以下の通り
1. build.pyを実行 (足りない残基を補う)
2. boxsize指定のため、system.pdb内を以下のように書き換える (100.000が3つ並んでいるのがboxsize)
```
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1 
```
3. run.shを実行 (sim.pyをsbatchで実行するスクリプト)