from setuptools import setup

setup(
    name='your-package',
    version='1.0',
    description='Your package description',
    install_requires=[
        'iqtree>=X.X.X',  # Specify the desired version of iqtree
        'qiime2-2022.2',
        'biom-format'
    ],
    packages=['your_package'],
)
