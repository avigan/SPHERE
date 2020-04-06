from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# setup
setup(
    name='vlt-sphere',
    version='1.1',
    description='Reduction and analysis code for the VLT/SPHERE instrument',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/avigan/SPHERE',
    author='Arthur Vigan',
    author_email='arthur.vigan@lam.fr',
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research ',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    keywords='vlt sphere ifs irdis reduction exoplanet pipeline',
    packages=['sphere', 'sphere.utils', 'sphere.IRDIS'],
    install_requires=[
        'numpy', 'scipy', 'astropy', 'pandas', 'matplotlib'
    ],
    include_package_data=True,
    package_data={
        'sphere': ['data/*.txt', 'data/*.dat', 'data/*.fits',
                   'instruments/*.ini', 'instruments/*.dat'],
    },
    zip_safe=False
)
