# Inputs for CHARMM36m alaine-dipeptide in explicit water molecules (TIP3P molecules)

Equilibration (NPTでの平衡化、 1 ns)

```shell
python 2_equilibration.py
```

Production (NPTでのプロダクション、 1 microsec)

```shell
python 3_production.py
```

Visualization with VMD

```shell
vmd -psf 1_setup/5_ionize/ionized.psf -dcd 3_production.dcd
```

Linux(GPU)で実行する場合と、Macで実行する場合、それぞれ `2_equilibration.py` と `3_production.py` の中身のどちらかのコメントを外してください（デフォルトでは以下のようにLinux(GPU)を使うようになっています)

```python
# Linux with GPU
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}

# Mac
#platform = Platform.getPlatformByName('OpenCL')
#platformProperties = {'Precision': 'single'}
```

バッチジョブを投入すると、 Equilibration と Productionがまとめて順番にバッチで実行されます。以下でバッチジョブを投入します。

```shell
sbatch --gres=gpu:a6000:1 -w floyd run.sh
```

