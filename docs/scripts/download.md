# Downloading data

Scripts for getting NICER data from the archive.

(get-nicer-aws)=
## get_nicer_aws.py

Finds the observations of a source with astroquery and downloads them from the HEASARC copy on Amazon Web Services. Needs the AWS command-line tool (`aws`).

```{argparse}
:filename: ../scripts/get_nicer_aws.py
:func: get_parser
:prog: get_nicer_aws.py
```

(ni-data-download)=
## ni_data_download.py

```{argparse}
:filename: ../scripts/ni_data_download.py
:func: get_parser
:prog: ni_data_download.py
```
