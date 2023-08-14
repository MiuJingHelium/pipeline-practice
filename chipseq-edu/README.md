# Into

A simple pipeline for educational purposes. The pipeline:
* aligns Chip-Seq single-end data
* generates reads coverage track

Peak calling step is missing on purposed and is supposed to be a simple exercise.

# Configure Pipeline
```shell
git clone https://github.com/JetBrains-Research/chipseq-smk-pipeline-edu#configure-pipeline <pipeline_working_dir>

cd <pipeline_working_dir>
```

Pipeline is tested with snakemake 6.3.0. It is recommended to create a fresh conda environment using `mamba` or `conda` (see https://snakemake.readthedocs.io/en/stable/getting_started/installation.html?highlight=mamba#installation-via-conda-mamba)
:
```shell
mamba env create --name snakemake --file ./environment.yaml
# or:
# conda env create --name snakemake --file ./environment.yaml
```
Activate conda environment with snakemake:
```shell
conda activate snakemake
```


**Configure input data**:
Download and unpack [example data](https://drive.google.com/file/d/1vtHVs4Yvf6ZnfynllxGYD7Lc96HuU46m/view?usp=sharing).
Copy `reads` folder to working directory. Delete `data_table.tsv` because `config` folder already has more recent version of `data_table.csv`.

# Run

Plot DAG and rule graphs
```shell
snakemake --dag | dot -Tsvg > images/dag.svg
snakemake --rulegraph | dot -Tsvg > images/rulegraph.svg
```

Check pipeline
```shell
snakemake -pr --use-conda --cores 8 --dry-run
```

Run pipeline
```shell
snakemake -pr --use-conda --cores 8

# or:
#snakemake -pr --use-conda --cores 8 --notemp
```

Archive pipeline results to a bundle
```shell
snakemake -pr --use-conda --cores 1 all_results_bundle
```

Some clusters automatically clean files older than 1 month. You could ask snakemake
touch all files in a correct pipeline-specific order. At the moment it doesn't work
with outputs marked as `temp(..)`, see https://github.com/snakemake/snakemake/issues/1028
```shell
./workflow/touch_pipeline.sh | grep -v "macs2" | grep -v ".idea"
```

##  Unit tests

### Run
```shell
pytest .tests/unit/test_bam_bigwig.py # Passes
pytest .tests/unit//test_reads_multiqc.py # Fails, custom matcher required
```
### Generate
Doesn't work smooth in Snakemake `6.5.1`, but could be done.

> **NB** Snakemake 6.5.1 tests generator doesn't work with 'multiqc' results (e.g. some case with directories in output), it fails with the error that directory already exists, workaround. Fix add `dirs_exist_ok=True`:
> ```python
>  # snakemake 6.5.1 : python3.8/site-packages/snakemake/unit_tests/__init__.py
>  # ~ line 81
>  if f.is_dir():
>      shutil.copytree(f, target, dirs_exist_ok=True)
>  ```


Assume the pipeline already executed pipeline:

```shell
 snakemake -pr --use-conda   --generate-unit-tests
```

Create config symlink:
```shell
cd .tests/unit
ln -s ../../config
```

Fix every text, add:
```python
config_path = PurePosixPath(".tests/unit/config")

# Copy data to the temporary workdir.
shutil.copytree(data_path, workdir)
shutil.copytree(config_path, workdir / "config") # copy config
```

Update `common.py`:
```python
# ...
for f in files:
    f = (Path(path) / f).relative_to(self.workdir)
    if (
            str(f).startswith(".snakemake")
            or str(f).startswith("config/") # add config
            or str(f).startswith("logs/")  # add logs
    ):
        continue
# ...
```