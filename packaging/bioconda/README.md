# Bioconda packaging

`meta.yaml` is a draft of the [Bioconda](https://bioconda.github.io/) recipe for MGCalibrator. It lives in this
repository for reference; Bioconda itself builds from the recipe in the
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes) repository.

## Release and submit

1. Make sure the version in `pyproject.toml` (and `.zenodo.json`) is the version to release, merge to `main`
   and create a GitHub release with the tag `v<version>` (e.g. `v1.1.0`).
2. Calculate the checksum of the release tarball:

   ```bash
   curl -sL https://github.com/NimroddeWit/MGCalibrator/archive/refs/tags/v1.1.0.tar.gz | sha256sum
   ```

3. Fork `bioconda/bioconda-recipes`, copy `meta.yaml` to `recipes/mgcalibrator/meta.yaml`, fill in the
   `sha256` and update `version` if needed.
4. Test the recipe locally (optional but recommended):

   ```bash
   mamba create -n bioconda-build -c conda-forge -c bioconda bioconda-utils
   mamba activate bioconda-build
   bioconda-utils lint --packages mgcalibrator recipes config.yml
   bioconda-utils build --packages mgcalibrator recipes config.yml
   ```

5. Open a pull request to `bioconda/bioconda-recipes`. After it is merged, the package is built and
   uploaded automatically and can be installed with:

   ```bash
   mamba install -c conda-forge -c bioconda mgcalibrator
   ```

## For new releases

After the first release, the Bioconda autobump bot opens pull requests for new GitHub releases. Keep the
dependencies in `pyproject.toml` and in the recipe in sync (and bump `build: number` for recipe-only changes).
