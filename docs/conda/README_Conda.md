Conda Environments Guide
=================

AVNI can be installed with various combinations of packages. We have tested the following
sets of packages. Please feel free to update these files or create new ones that work
for specific applications.
* `environment.yml` — Base version that all public CIG users should use. Excludes modules
for documentation and applets.
* `environment_doc.yml` — The version for creating documentation and testing through
Travis CI includes doc modules such as sphinx but exclude applets.
* `environment_devel.yml` — The version that all core developers and maintainers should use. Most comprehensive and up to date.
* `environment_atlas3d.yml` — The version for model assimilation on tiger/tigressdata
