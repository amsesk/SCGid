from setuptools import setup
setup(
    name = "SCGid",
    version = "0.9b",
    author = "Kevin Amses",
    author_email = "amsesk@umich.edu",
    description = (
        "A consensus approach to contig filtering and genome prediction from single-cell sequencing libraries"
        ),
    license = "GPL-3.0",
    keywords = "",
    url = "http://www.github.com/amsesk/SCGid",
    packages=["scgid"],
    data_files=[
        ('scgid', ['scgid/print_tetramer_freqs_deg_filter_esom_VD.pl', 'scgid/ape_nj.R', 'scgid/codons_phytools.R', 'scgid/gc_cov.plot.R']),
    ],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'scgid = scgid.main:main'
        ]
    },
    install_requires=[
        'numpy>=1.15.0',
        'pandas>=0.23.4',
        'ete3>=3.1.1',
        'plotly',
        'pyyaml',
        'urllib3',
        'importlib_metadata'
    ]
)
