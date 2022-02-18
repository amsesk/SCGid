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
    entry_points={
        'console_scripts': [
            'scgid = scgid.main:main'
        ]
    },
    setup_requires=[
        'numpy>=1.15.0',
        'cython==0.29.28'
        ],
    install_requires=[
        'numpy>=1.15.0',
        'cython==0.29.28',
        'pandas>=0.23.4',
        'ete3>=3.1.1',
        'plotly',
        'pyyaml',
        'urllib3',
        'importlib_metadata'
    ]
)
