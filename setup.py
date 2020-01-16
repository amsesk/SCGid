from setuptools import setup
setup(
    name = "SCGid",
    version = "1.0.0",
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
        ('scgid', ['scgid/logging_config.ini', 'scgid/print_tetramer_freqs_deg_filter_esom_VD.pl']),
    ],
    include_package_data=True,
    #py_modules=['scgid.scgid', 'scgid.gct', 'scgid.codons', 'scgid.kmers', 'logging'],
    #scripts=['bin/scgid'],
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

# Import config script after above
